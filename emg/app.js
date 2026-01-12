/**
 * EMG walking app (serverless, browser-only) with optional VIDEO assistance.
 *
 * Key policy:
 *  - CSV1 is TA (recommended). HS is estimated from TA.
 *  - CSV2 is the muscle of interest. Its gait cycle (0-100%) is segmented using TA HS.
 *  - Video (mp4) can help select the correct HS sequence when TA has double bursts per cycle.
 *
 * Processing:
 *  - Bandpass: HP then LP (2nd-order biquad applied twice, filtfilt-like)
 *  - Rectify
 *  - RMS envelope (centered)
 *  - HS candidates: burst onset by MAD threshold + peak fallback
 *  - TA HS refinement:
 *      - remove events before baselineSec (e.g., 3s standing)
 *      - if video assist ON and video events exist: auto-sync + keep matched HS
 *      - else: every-other decimation if looks doubled
 *  - %peak normalization
 *  - Cycle normalization to N points (0-100%)
 *  - Plot mean±SD (+ individual)
 */

const $ = (id) => document.getElementById(id);

const state = {
  file1: null,
  file2: null,
  parsed1: null,
  parsed2: null,
  result1: null,
  result2: null,

  // video
  videoFile: null,
  roi: null,            // {cx, cy, size} in VIDEO SOURCE coordinates
  videoEvents: null,    // [sec]
  videoEnergy: null,    // {t:[], e:[], thr:number}
};

function setMsg(s, isErr = false) {
  const el = $("msg");
  el.textContent = s;
  el.style.color = isErr ? "#ffb4b4" : "rgba(233,238,255,.72)";
}
function setVideoStatus(s, isErr = false) {
  const el = $("videoStatus");
  el.textContent = s;
  el.style.color = isErr ? "#ffb4b4" : "rgba(233,238,255,.72)";
}

function median(arr) {
  const a = Array.from(arr).filter(Number.isFinite);
  if (a.length === 0) return NaN;
  a.sort((x, y) => x - y);
  const mid = Math.floor(a.length / 2);
  return a.length % 2 === 0 ? (a[mid - 1] + a[mid]) / 2 : a[mid];
}
function mad(arr) {
  const med = median(arr);
  const dev = arr.map((v) => Math.abs(v - med));
  return median(dev);
}
function percentile(arr, p) {
  const a = Array.from(arr).filter(Number.isFinite).sort((x, y) => x - y);
  if (a.length === 0) return NaN;
  const pos = (p / 100) * (a.length - 1);
  const i0 = Math.floor(pos);
  const i1 = Math.min(a.length - 1, i0 + 1);
  const w = pos - i0;
  return a[i0] * (1 - w) + a[i1] * w;
}
function linspace(a, b, n) {
  const out = new Float64Array(n);
  const step = (b - a) / (n - 1);
  for (let i = 0; i < n; i++) out[i] = a + step * i;
  return out;
}

function reverseArray(x) {
  const y = new Float64Array(x.length);
  for (let i = 0; i < x.length; i++) y[i] = x[x.length - 1 - i];
  return y;
}

function demean(x) {
  let s = 0, n = 0;
  for (const v of x) if (Number.isFinite(v)) { s += v; n++; }
  const m = n > 0 ? s / n : 0;
  const y = new Float64Array(x.length);
  for (let i = 0; i < x.length; i++) y[i] = x[i] - m;
  return y;
}

function absArray(x) {
  const y = new Float64Array(x.length);
  for (let i = 0; i < x.length; i++) y[i] = Math.abs(x[i]);
  return y;
}

/** centered moving RMS */
function movingRMS(x, win) {
  const n = x.length;
  const half = Math.floor(win / 2);
  const sq = new Float64Array(n);
  for (let i = 0; i < n; i++) sq[i] = x[i] * x[i];

  const prefix = new Float64Array(n + 1);
  for (let i = 0; i < n; i++) prefix[i + 1] = prefix[i] + sq[i];

  const out = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    const s = Math.max(0, i - half);
    const e = Math.min(n, i + half + 1);
    const mean = (prefix[e] - prefix[s]) / (e - s);
    out[i] = Math.sqrt(mean);
  }
  return out;
}

/** RBJ biquad coefficients (Butterworth Q=1/sqrt(2)) */
function biquadCoeffs(type, fc, fs, Q = 1 / Math.SQRT2) {
  const w0 = (2 * Math.PI * fc) / fs;
  const cosw0 = Math.cos(w0);
  const sinw0 = Math.sin(w0);
  const alpha = sinw0 / (2 * Q);

  let b0, b1, b2, a0, a1, a2;
  if (type === "lowpass") {
    b0 = (1 - cosw0) / 2;
    b1 = 1 - cosw0;
    b2 = (1 - cosw0) / 2;
    a0 = 1 + alpha;
    a1 = -2 * cosw0;
    a2 = 1 - alpha;
  } else if (type === "highpass") {
    b0 = (1 + cosw0) / 2;
    b1 = -(1 + cosw0);
    b2 = (1 + cosw0) / 2;
    a0 = 1 + alpha;
    a1 = -2 * cosw0;
    a2 = 1 - alpha;
  } else {
    throw new Error("unknown biquad type");
  }
  return {
    b0: b0 / a0,
    b1: b1 / a0,
    b2: b2 / a0,
    a1: a1 / a0,
    a2: a2 / a0,
  };
}

function biquadFilter(x, c) {
  const y = new Float64Array(x.length);
  let x1 = 0, x2 = 0, y1 = 0, y2 = 0;
  for (let i = 0; i < x.length; i++) {
    const x0 = x[i];
    const y0 = c.b0 * x0 + c.b1 * x1 + c.b2 * x2 - c.a1 * y1 - c.a2 * y2;
    y[i] = y0;
    x2 = x1; x1 = x0;
    y2 = y1; y1 = y0;
  }
  return y;
}

/** forward-backward filtering (approx) */
function filtfilt(x, coeffs) {
  let y = biquadFilter(x, coeffs);
  y = reverseArray(y);
  y = biquadFilter(y, coeffs);
  y = reverseArray(y);
  return y;
}

/** 4th-order-ish by applying same biquad twice */
function butterworth2xFILT(x, type, fc, fs) {
  const c = biquadCoeffs(type, fc, fs, 1 / Math.SQRT2);
  let y = filtfilt(x, c);
  y = filtfilt(y, c);
  return y;
}

function risingEdges(mask) {
  const starts = [];
  const ends = [];
  for (let i = 1; i < mask.length; i++) {
    if (!mask[i - 1] && mask[i]) starts.push(i);
    if (mask[i - 1] && !mask[i]) ends.push(i);
  }
  if (mask[0]) starts.unshift(0);
  if (mask[mask.length - 1]) ends.push(mask.length);

  const segs = [];
  for (let i = 0; i < Math.min(starts.length, ends.length); i++) segs.push([starts[i], ends[i]]);
  return segs;
}

function localMaxima(y) {
  const peaks = [];
  for (let i = 1; i < y.length - 1; i++) {
    if (y[i] > y[i - 1] && y[i] >= y[i + 1]) peaks.push(i);
  }
  return peaks;
}

/** Peak picking with min distance */
function pickPeaksGreedy(y, minDist, thr, maxCount = null) {
  const peaks = localMaxima(y).filter((i) => y[i] > thr);
  peaks.sort((i, j) => y[j] - y[i]); // amplitude desc
  const chosen = [];
  for (const idx of peaks) {
    let ok = true;
    for (const c of chosen) {
      if (Math.abs(idx - c) < minDist) { ok = false; break; }
    }
    if (ok) {
      chosen.push(idx);
      if (maxCount && chosen.length >= maxCount) break;
    }
  }
  chosen.sort((a, b) => a - b);
  return chosen;
}

function refineToOnset(y, peakIdx, thrOff) {
  let i = peakIdx;
  while (i > 0 && y[i] > thrOff) i--;
  return i;
}

function scoreHS(hs, fs, plausibleStep, expectedCount=null) {
  if (!hs || hs.length < 2) return Infinity;
  const intervals = [];
  for (let i = 1; i < hs.length; i++) intervals.push((hs[i] - hs[i - 1]) / fs);

  const meanInt = intervals.reduce((a, b) => a + b, 0) / intervals.length;
  const medInt = median(intervals);
  const sdInt = Math.sqrt(intervals.reduce((a, b) => a + (b - meanInt) * (b - meanInt), 0) / intervals.length);
  const cv = sdInt / (meanInt + 1e-12);
  const [lo, hi] = plausibleStep;
  const bad = intervals.filter((v) => v < lo || v > hi).length / intervals.length;

  let score = 5.0 * bad + 2.0 * cv + 0.5 * Math.abs(medInt - 0.7);
  if (expectedCount !== null) score += 0.8 * Math.abs(hs.length - expectedCount);
  return score;
}

function applyOffsetAndMerge(hsRaw, offset, minGap, nMax) {
  if (!hsRaw || hsRaw.length === 0) return [];
  let hs = hsRaw
    .map((v) => v + offset)
    .filter((v) => v >= 0 && v < nMax)
    .sort((a, b) => a - b);

  const merged = [hs[0]];
  for (let i = 1; i < hs.length; i++) {
    if (hs[i] - merged[merged.length - 1] >= minGap) merged.push(hs[i]);
  }
  return merged;
}

/**
 * HS detection from envelope:
 *  - primary: MAD threshold burst-onset
 *  - fallback: peak-based
 */
function detectHS(env, fs, opts) {
  const {
    expectedCount = null,
    minBurstMs = 30,
    minGapMs = 450,
    offsetMs = 80,
    plausibleStep = [0.35, 1.6],
  } = opts;

  const minBurst = Math.round((fs * minBurstMs) / 1000);
  const minGap = Math.round((fs * minGapMs) / 1000);
  const offset = Math.round((fs * offsetMs) / 1000);

  const med = median(env);
  let m = mad(env);

  if (!Number.isFinite(m) || m < 1e-12) {
    const p75 = percentile(env, 75);
    const p25 = percentile(env, 25);
    m = (p75 - p25) / 1.349;
  }

  const ks = [];
  for (let k = 1.2; k <= 4.0 + 1e-9; k += 0.1) ks.push(Math.round(k * 100) / 100);

  let best = null;

  // onset-based
  if (Number.isFinite(m) && m > 1e-12) {
    for (const k of ks) {
      const thr = med + k * m;
      const mask = new Array(env.length);
      for (let i = 0; i < env.length; i++) mask[i] = env[i] > thr;

      let segs = risingEdges(mask);
      segs = segs.filter(([s, e]) => (e - s) >= minBurst);
      if (segs.length === 0) continue;

      let hs = segs.map(([s, _e]) => s);
      hs = applyOffsetAndMerge(hs, offset, minGap, env.length);
      if (hs.length < 2) continue;

      const sc = scoreHS(hs, fs, plausibleStep, expectedCount);
      if (best === null || sc < best.score) best = { k, thr, hs, score: sc };
    }
  }
  if (best) return best;

  // peak fallback
  const thrPeak = percentile(env, 80);
  const thrOff = percentile(env, 60);
  let peaks = pickPeaksGreedy(env, minGap, thrPeak, expectedCount);
  if (peaks.length >= 2) {
    const onsets = peaks.map((p) => refineToOnset(env, p, thrOff));
    const hs = applyOffsetAndMerge(onsets, offset, minGap, env.length);
    if (hs.length >= 2) return { k: "peak-fallback", thr: thrPeak, hs, score: scoreHS(hs, fs, plausibleStep, expectedCount) };
  }

  peaks = applyOffsetAndMerge(peaks, offset, minGap, env.length);
  if (peaks.length >= 2) return { k: "peak-direct", thr: thrPeak, hs: peaks, score: scoreHS(peaks, fs, plausibleStep, expectedCount) };

  throw new Error("HS推定に失敗しました（信号が弱い/ノイズ多い/閾値設定不適）。minGapやHP(20-30Hz)やオフセットを調整してください。");
}

function percentPeak(env, hs) {
  let ref = NaN;
  if (hs && hs.length >= 2) {
    const s = hs[0], e = hs[hs.length - 1];
    let mx = -Infinity;
    for (let i = s; i < e; i++) if (env[i] > mx) mx = env[i];
    ref = mx;
  } else {
    ref = Math.max(...env);
  }
  if (!Number.isFinite(ref) || ref <= 0) return new Float64Array(env.length).fill(NaN);

  const out = new Float64Array(env.length);
  for (let i = 0; i < env.length; i++) out[i] = (env[i] / ref) * 100;
  return out;
}

function normalizeCycles(y, hs, nPoints) {
  const grid = linspace(0, 100, nPoints);
  const cycles = [];
  for (let i = 0; i < hs.length - 1; i++) {
    const s = hs[i], e = hs[i + 1];
    if (e <= s + 2) continue;

    const segLen = e - s;
    const yseg = y.subarray(s, e);
    const out = new Float64Array(nPoints);
    const last = segLen - 1;

    for (let j = 0; j < nPoints; j++) {
      const pos = (j / (nPoints - 1)) * last;
      const i0 = Math.floor(pos);
      const i1 = Math.min(last, i0 + 1);
      const w = pos - i0;
      out[j] = yseg[i0] * (1 - w) + yseg[i1] * w;
    }
    cycles.push(out);
  }
  return { grid, cycles };
}

function meanAndSD(cycles) {
  const n = cycles.length;
  if (n === 0) return { mean: [], sd: [] };

  const m = cycles[0].length;
  const mean = new Float64Array(m);
  const sd = new Float64Array(m);

  for (let j = 0; j < m; j++) {
    let s = 0, c = 0;
    for (let i = 0; i < n; i++) {
      const v = cycles[i][j];
      if (Number.isFinite(v)) { s += v; c++; }
    }
    mean[j] = c > 0 ? s / c : NaN;
  }

  for (let j = 0; j < m; j++) {
    let s = 0, c = 0;
    for (let i = 0; i < n; i++) {
      const v = cycles[i][j];
      if (Number.isFinite(v)) { s += (v - mean[j]) * (v - mean[j]); c++; }
    }
    sd[j] = c > 1 ? Math.sqrt(s / c) : NaN;
  }
  return { mean, sd };
}

/* =======================
   iPhone Safari対策：CSV読み取り
   ======================= */

async function readFileData(file) {
  const buf = await file.arrayBuffer();
  let text = "";
  try {
    text = new TextDecoder("utf-8", { fatal: false }).decode(buf);
  } catch (e) {
    text = "";
  }
  return { buf, text };
}
function bytesToLatin1String(buf) {
  try {
    return new TextDecoder("iso-8859-1", { fatal: false }).decode(buf);
  } catch (e) {
    const bytes = new Uint8Array(buf);
    let out = "";
    const CHUNK = 4096;
    for (let i = 0; i < bytes.length; i += CHUNK) {
      const sub = bytes.subarray(i, i + CHUNK);
      out += String.fromCharCode.apply(null, sub);
    }
    return out;
  }
}

function parseCSV(buf, text, filename, fsDefault) {
  const name = filename.replace(/\.csv$/i, "");

  const raw = bytesToLatin1String(buf);
  const rawLines = raw.split(/\r\n|\n|\r/);

  // sensor format: find "No,"
  let start = -1;
  let headerLine = "";
  for (let i = 0; i < Math.min(rawLines.length, 5000); i++) {
    const line = rawLines[i];
    if (!line) continue;
    const idx = line.indexOf("No,");
    if (idx >= 0) {
      start = i;
      headerLine = line.slice(idx);
      break;
    }
  }

  if (start >= 0) {
    const header = headerLine.split(",");
    const idxNo = header.findIndex((h) => (h || "").trim() === "No");

    let idxV = header.findIndex((h) =>
      (h || "").includes("mV") ||
      (h || "").toLowerCase().includes("voltage") ||
      (h || "").includes("電圧")
    );
    if (idxV < 0) idxV = 1;

    const no = [];
    const x = [];

    for (let i = start + 1; i < rawLines.length; i++) {
      const line = rawLines[i];
      if (!line) continue;
      const parts = line.split(",");
      if (parts.length < Math.max(idxNo, idxV) + 1) continue;

      const n = parseInt(parts[idxNo], 10);
      const v = parseFloat(parts[idxV]);
      if (Number.isFinite(n) && Number.isFinite(v)) {
        no.push(n);
        x.push(v);
      }
    }

    const t = new Float64Array(no.length);
    for (let i = 0; i < no.length; i++) t[i] = (no[i] - 1) / fsDefault;

    return {
      name,
      fs: fsDefault,
      t,
      x: Float64Array.from(x),
      columns: ["No", "value"],
      guessedTimeCol: "No",
      guessedSigCol: "value",
      table: null,
      format: "sensor",
    };
  }

  // table format
  let lines = (text || "").split(/\r\n|\n|\r/);
  if (!lines[0] || lines[0].split(",").length < 2) {
    lines = rawLines.filter((l) => l && l.includes(","));
  }

  const header = (lines[0] || "").split(",").map((s) => (s || "").trim()).filter(Boolean);

  const rows = [];
  for (let i = 1; i < lines.length; i++) {
    const line = lines[i];
    if (!line) continue;
    const parts = line.split(",");
    if (parts.length < header.length) continue;
    rows.push(parts.slice(0, header.length));
  }

  const cols = {};
  for (let j = 0; j < header.length; j++) {
    const key = header[j];
    cols[key] = new Float64Array(rows.length);
    for (let i = 0; i < rows.length; i++) cols[key][i] = parseFloat(rows[i][j]);
  }

  let timeCol = header.find((h) => ["time", "time_s", "t", "sec", "seconds"].includes((h || "").toLowerCase())) || null;
  if (!timeCol) timeCol = header[0];

  let sigCol = header.find((h) => h !== timeCol) || header[0];

  return {
    name,
    fs: fsDefault,
    t: cols[timeCol],
    x: cols[sigCol],
    columns: header,
    guessedTimeCol: timeCol,
    guessedSigCol: sigCol,
    table: cols,
    format: "table",
  };
}

function populateSelect(sel, options, selected) {
  sel.innerHTML = "";
  for (const opt of options) {
    const o = document.createElement("option");
    o.value = opt;
    o.textContent = opt;
    if (opt === selected) o.selected = true;
    sel.appendChild(o);
  }
}

function getSelectedSignal(parsed, sigSel, timeSel, fs) {
  if (parsed.format === "sensor") return { t: parsed.t, x: parsed.x };

  const sig = sigSel.value || parsed.guessedSigCol;
  const tim = timeSel.value || parsed.guessedTimeCol;

  const x = parsed.table[sig];
  const t0 = parsed.table[tim];

  const maxT = Math.max(...t0);
  const estDur = t0.length / fs;

  // if time column looks like seconds, keep
  const isSeconds = Number.isFinite(maxT) && maxT > estDur * 0.8;
  if (isSeconds) return { t: t0, x };

  const tSec = new Float64Array(t0.length);
  for (let i = 0; i < t0.length; i++) tSec[i] = (t0[i] - t0[0]) / fs;
  return { t: tSec, x };
}

/* =======================
   Compute pipeline
   ======================= */

function computePipeline(t, x, fs, params, expectedCount) {
  const { hp, lp, rmsMs, minBurstMs, minGapMs, offsetMs, nPoints } = params;

  let y = demean(x);

  if (hp > 0) y = butterworth2xFILT(y, "highpass", hp, fs);
  if (lp > 0 && lp < fs / 2) y = butterworth2xFILT(y, "lowpass", lp, fs);

  y = absArray(y);
  const win = Math.max(1, Math.round((fs * rmsMs) / 1000));
  const env = movingRMS(y, win);

  const hsRes = detectHS(env, fs, {
    expectedCount,
    minBurstMs,
    minGapMs,
    offsetMs,
  });

  const envPct = percentPeak(env, hsRes.hs);
  const cyc = normalizeCycles(envPct, hsRes.hs, nPoints);
  const ms = meanAndSD(cyc.cycles);

  return {
    t,
    envRaw: env,
    envPct,
    hs: hsRes.hs,
    k: hsRes.k,
    thr: hsRes.thr,
    grid: cyc.grid,
    cycles: cyc.cycles,
    mean: ms.mean,
    sd: ms.sd,
  };
}

function applyHsToExistingResult(res, hs, nPoints) {
  const hs2 = hs.filter((i) => i >= 0 && i < res.envRaw.length);
  const envPct = percentPeak(res.envRaw, hs2);
  const cyc = normalizeCycles(envPct, hs2, nPoints);
  const ms = meanAndSD(cyc.cycles);
  return { ...res, hs: hs2, envPct, grid: cyc.grid, cycles: cyc.cycles, mean: ms.mean, sd: ms.sd };
}

// if TA looks doubled, choose odd/even by score
function chooseEveryOtherIfDouble(hs, fs, expectedCount = null) {
  if (!hs || hs.length < 6) return hs;

  const intervals = [];
  for (let i = 1; i < hs.length; i++) intervals.push((hs[i] - hs[i - 1]) / fs);
  const medInt = median(intervals);

  const looksDouble = (medInt < 0.50) || (expectedCount !== null && hs.length >= expectedCount * 1.6);
  if (!looksDouble) return hs;

  const candA = hs.filter((_, i) => i % 2 === 0);
  const candB = hs.filter((_, i) => i % 2 === 1);

  const sA = scoreHS(candA, fs, [0.35, 1.6], expectedCount);
  const sB = scoreHS(candB, fs, [0.35, 1.6], expectedCount);
  return (sA <= sB) ? candA : candB;
}

/* =======================
   Plotting
   ======================= */

function plotCycle(targetId, res, showIndividual) {
  const traces = [];

  if (showIndividual) {
    for (const c of res.cycles) {
      traces.push({
        x: Array.from(res.grid),
        y: Array.from(c),
        mode: "lines",
        line: { width: 1, color: "rgba(200,200,200,0.35)" },
        hoverinfo: "skip",
        showlegend: false,
      });
    }
  }

  if (res.cycles.length > 0) {
    const x = Array.from(res.grid);
    const mean = Array.from(res.mean);
    const sd = Array.from(res.sd);
    const upper = mean.map((v, i) => v + sd[i]);
    const lower = mean.map((v, i) => v - sd[i]);

    traces.push({ x, y: upper, mode: "lines", line: { width: 0 }, hoverinfo: "skip", showlegend: false });
    traces.push({
      x,
      y: lower,
      mode: "lines",
      fill: "tonexty",
      fillcolor: "rgba(59,130,246,0.18)",
      line: { width: 0 },
      hoverinfo: "skip",
      showlegend: false,
    });
    traces.push({
      x,
      y: mean,
      mode: "lines",
      line: { width: 3, color: "#3b82f6" },
      name: "mean",
    });
  }

  const layout = {
    margin: { l: 55, r: 20, t: 30, b: 45 },
    xaxis: { title: "gait cycle (%)", range: [0, 100] },
    yaxis: { title: "%RMS (%peak)" },
    paper_bgcolor: "rgba(0,0,0,0)",
    plot_bgcolor: "rgba(0,0,0,0)",
    font: { color: "#e9eeff" },
  };

  return Plotly.newPlot(targetId, traces, layout, { displayModeBar: false, responsive: true });
}

function plotDebug(targetId, res) {
  const t = Array.from(res.t);
  const y = Array.from(res.envPct);

  const shapes = (res.hs || []).map((idx) => ({
    type: "line",
    x0: res.t[idx],
    x1: res.t[idx],
    y0: 0,
    y1: 1,
    xref: "x",
    yref: "paper",
    line: { width: 1, color: "rgba(255,255,255,0.25)" },
  }));

  const traces = [{
    x: t,
    y: y,
    mode: "lines",
    line: { width: 2, color: "#22c55e" },
    name: "%peak env",
  }];

  const layout = {
    margin: { l: 55, r: 20, t: 30, b: 45 },
    xaxis: { title: "time (s)" },
    yaxis: { title: "%RMS (%peak)" },
    shapes,
    paper_bgcolor: "rgba(0,0,0,0)",
    plot_bgcolor: "rgba(0,0,0,0)",
    font: { color: "#e9eeff" },
  };

  return Plotly.newPlot(targetId, traces, layout, { displayModeBar: false, responsive: true });
}

function plotVideoDebug(targetId, tArr, eArr, thr, events) {
  const shapes = (events || []).map((ts) => ({
    type: "line",
    x0: ts, x1: ts,
    y0: 0, y1: 1,
    xref: "x", yref: "paper",
    line: { width: 1, color: "rgba(255,255,255,0.25)" },
  }));

  const traces = [{
    x: tArr,
    y: eArr,
    mode: "lines",
    line: { width: 2, color: "#f59e0b" },
    name: "motion energy",
  }];

  if (Number.isFinite(thr)) {
    traces.push({
      x: [tArr[0], tArr[tArr.length-1]],
      y: [thr, thr],
      mode: "lines",
      line: { width: 2, color: "rgba(59,130,246,0.7)", dash: "dash" },
      name: "thr",
    });
  }

  const layout = {
    margin: { l: 55, r: 20, t: 30, b: 45 },
    xaxis: { title: "video time (s)" },
    yaxis: { title: "motion energy (a.u.)" },
    shapes,
    paper_bgcolor: "rgba(0,0,0,0)",
    plot_bgcolor: "rgba(0,0,0,0)",
    font: { color: "#e9eeff" },
  };
  return Plotly.newPlot(targetId, traces, layout, { displayModeBar: false, responsive: true });
}

/* =======================
   Video assist
   ======================= */
function fileToDataURL(file){
  return new Promise((resolve, reject) => {
    const r = new FileReader();
    r.onload = () => resolve(r.result);
    r.onerror = () => reject(r.error || new Error("FileReader error"));
    r.readAsDataURL(file);
  });
}

function waitLoadedMetadata(video, timeoutMs){
  return new Promise((resolve) => {
    let done = false;
    const timer = setTimeout(() => {
      if (done) return;
      done = true;
      cleanup();
      resolve(false);
    }, timeoutMs);

    function cleanup(){
      clearTimeout(timer);
      video.removeEventListener("loadedmetadata", onOk);
      video.removeEventListener("error", onErr);
    }
    function onOk(){
      if (done) return;
      done = true;
      cleanup();
      resolve(true);
    }
    function onErr(){
      if (done) return;
      done = true;
      cleanup();
      resolve(false);
    }

    video.addEventListener("loadedmetadata", onOk, { once:true });
    video.addEventListener("error", onErr, { once:true });
  });
}

function videoErrorText(video){
  const e = video?.error;
  if(!e) return "unknown";
  // 1:ABORTED 2:NETWORK 3:DECODE 4:SRC_NOT_SUPPORTED
  const map = {1:"ABORTED",2:"NETWORK",3:"DECODE",4:"SRC_NOT_SUPPORTED"};
  return `${map[e.code] || e.code}`;
}

async function loadVideoFileToElement(file){
  const video = $("video");

  // cleanup old
  try { video.pause(); } catch {}
  try {
    if (video.src && video.src.startsWith("blob:")) URL.revokeObjectURL(video.src);
  } catch {}
  video.removeAttribute("src");
  video.load();

  // ensure inline
  video.muted = true;
  video.playsInline = true;
  video.setAttribute("playsinline", "");
  video.setAttribute("webkit-playsinline", "");
  video.preload = "metadata";

  // 1) try blob URL (fast)
  const blobUrl = URL.createObjectURL(file);
  video.src = blobUrl;
  video.load();

  let ok = await waitLoadedMetadata(video, 1500);
  if(!ok){
    // 2) fallback to DataURL (more compatible on iOS, but memory heavier)
    try { URL.revokeObjectURL(blobUrl); } catch {}
    const dataUrl = await fileToDataURL(file);
    video.removeAttribute("src");
    video.load();
    video.src = dataUrl;
    video.load();

    ok = await waitLoadedMetadata(video, 5000);
    if(!ok){
      throw new Error(`動画を読み込めません（${videoErrorText(video)}）。iPhoneでは高fps動画が失敗することがあります。30/60fpsのH.264で書き出して再試行してください。`);
    }
  }

  // show first frame (sometimes black until seek)
  try{
    const t = Math.min(0.05, Math.max(0, (video.duration || 0) - 0.05));
    video.currentTime = t;
    video.pause();
  } catch {}
}

function drawRoiOverlay() {
  const video = $("video");
  const canvas = $("overlay");
  const ctx = canvas.getContext("2d");

  // match display size
  const rect = video.getBoundingClientRect();
  canvas.width = Math.max(1, Math.round(rect.width * devicePixelRatio));
  canvas.height = Math.max(1, Math.round(rect.height * devicePixelRatio));
  ctx.setTransform(devicePixelRatio, 0, 0, devicePixelRatio, 0, 0);

  ctx.clearRect(0, 0, rect.width, rect.height);

  if (!state.roi || !video.videoWidth || !video.videoHeight) return;

  const size = parseFloat($("roiSize").value);
  const scaleX = rect.width / video.videoWidth;
  const scaleY = rect.height / video.videoHeight;

  const x = (state.roi.cx - size / 2) * scaleX;
  const y = (state.roi.cy - size / 2) * scaleY;
  const w = size * scaleX;
  const h = size * scaleY;

  ctx.strokeStyle = "rgba(34,197,94,0.9)";
  ctx.lineWidth = 3;
  ctx.strokeRect(x, y, w, h);
  ctx.fillStyle = "rgba(34,197,94,0.12)";
  ctx.fillRect(x, y, w, h);
}

function setRoiFromTap(clientX, clientY) {
  const video = $("video");
  const rect = video.getBoundingClientRect();
  if (!video.videoWidth || !video.videoHeight) return;

  const x = (clientX - rect.left) / rect.width;  // 0..1
  const y = (clientY - rect.top) / rect.height;  // 0..1

  state.roi = {
    cx: x * video.videoWidth,
    cy: y * video.videoHeight,
  };
  drawRoiOverlay();
}

function smoothMovingAvg(y, win) {
  if (win <= 1) return y.slice();
  const out = new Float64Array(y.length);
  let s = 0;
  let q = [];
  for (let i = 0; i < y.length; i++) {
    s += y[i];
    q.push(y[i]);
    if (q.length > win) s -= q.shift();
    out[i] = s / q.length;
  }
  return out;
}

async function extractVideoEvents() {
  const video = $("video");
  if (!video.src) throw new Error("動画が未選択です。");
  if (!state.roi) throw new Error("動画を1回タップしてROI（足付近）を指定してください。");

  const roiSize = parseFloat($("roiSize").value);
  const baselineSec = parseFloat($("baselineSec").value);
  const minGapMs = parseFloat($("videoMinGap").value);
  const fpsTarget = Math.max(5, Math.min(30, parseFloat($("videoFps").value)));

  // analysis canvas (downsample ROI)
  const W = 64, H = 64;
  const c = document.createElement("canvas");
  c.width = W; c.height = H;
  const ctx = c.getContext("2d", { willReadFrequently: true });

  const prev = new Float32Array(W * H);
  let hasPrev = false;

  const times = [];
  const energy = [];

  // clamp ROI in source coordinates
  const size = Math.max(40, Math.min(Math.min(video.videoWidth, video.videoHeight), roiSize));
  const sx0 = Math.max(0, Math.min(video.videoWidth - size, state.roi.cx - size / 2));
  const sy0 = Math.max(0, Math.min(video.videoHeight - size, state.roi.cy - size / 2));

  // play & capture frames
  let lastT = -1;

  async function captureWithRVFC() {
    return new Promise((resolve, reject) => {
      let stopped = false;

      function stopResolve() {
        if (stopped) return;
        stopped = true;
        try { video.pause(); } catch {}
        resolve();
      }

      video.currentTime = 0;
      video.muted = true;
      video.playbackRate = 1.0;

      const stepMin = 1 / fpsTarget;

      const onFrame = (_now, meta) => {
        if (stopped) return;

        const t = meta.mediaTime;
        if (lastT < 0 || (t - lastT) >= stepMin) {
          lastT = t;

          ctx.drawImage(video, sx0, sy0, size, size, 0, 0, W, H);
          const img = ctx.getImageData(0, 0, W, H).data;

          let sum = 0;
          let idx = 0;
          // sample every pixel (still OK at 4096)
          for (let i = 0; i < img.length; i += 4) {
            const r = img[i], g = img[i + 1], b = img[i + 2];
            const lum = (0.299 * r + 0.587 * g + 0.114 * b) / 255.0;
            if (hasPrev) sum += Math.abs(lum - prev[idx]);
            prev[idx] = lum;
            idx++;
          }

          times.push(t);
          energy.push(hasPrev ? (sum / (W * H)) : 0);
          hasPrev = true;
        }

        if (t >= video.duration || video.ended) {
          stopResolve();
          return;
        }
        if (video.requestVideoFrameCallback) {
          video.requestVideoFrameCallback(onFrame);
        } else {
          stopResolve();
        }
      };

      video.play().then(() => {
        if (!video.requestVideoFrameCallback) {
          reject(new Error("このブラウザはrequestVideoFrameCallback非対応です。別ブラウザ/OSでお試しください。"));
          return;
        }
        video.requestVideoFrameCallback(onFrame);
      }).catch((e) => reject(e));

      video.onended = stopResolve;
    });
  }

  if (!video.requestVideoFrameCallback) {
    throw new Error("このブラウザは動画フレームコールバックに未対応です（iOSが古い可能性）。");
  }

  setVideoStatus("動画解析中…（自動再生します）");
  await captureWithRVFC();

  if (times.length < 10) throw new Error("動画フレームが十分取得できませんでした。");

  // smooth energy
  const dt = (times[times.length - 1] - times[0]) / Math.max(1, times.length - 1);
  const win = Math.max(1, Math.round(0.15 / Math.max(1e-3, dt)));
  const eSmooth = smoothMovingAvg(Float64Array.from(energy), win);

  // baseline threshold from first baselineSec seconds
  const base = [];
  for (let i = 0; i < times.length; i++) if (times[i] <= baselineSec) base.push(eSmooth[i]);
  if (base.length < 5) {
    for (let i = 0; i < Math.min(times.length, 60); i++) base.push(eSmooth[i]);
  }
  const medE = median(base);
  let mE = mad(base);
  if (!Number.isFinite(mE) || mE < 1e-9) mE = (percentile(base, 75) - percentile(base, 25)) / 1.349;
  const thr = medE + 4.0 * mE; // fixed k=4 for video motion

  const mask = new Array(times.length);
  for (let i = 0; i < times.length; i++) mask[i] = eSmooth[i] > thr;

  // event = end of "high motion" segment (falling edge) => foot stops (HS-ish)
  const events = [];
  const minGapSec = minGapMs / 1000;

  for (let i = 1; i < mask.length; i++) {
    if (mask[i - 1] && !mask[i]) {
      const tEvent = times[i]; // crossing time
      if (tEvent < baselineSec) continue;
      if (events.length === 0 || (tEvent - events[events.length - 1]) >= minGapSec) {
        events.push(tEvent);
      }
    }
  }

  state.videoEvents = events;
  state.videoEnergy = { t: times.slice(), e: Array.from(eSmooth), thr };

  plotVideoDebug("videoPlot", times, Array.from(eSmooth), thr, events);

  setVideoStatus(`動画イベント抽出OK：${events.length}個（ROI=緑枠）`);
  return events;
}

// find best time shift between EMG times and video events
function bestShiftMatch(emgTimes, videoTimes, rangeSec, stepSec, tolSec) {
  if (!emgTimes.length || !videoTimes.length) return { shift: 0, matchIdx: [], score: -1 };

  const vt = videoTimes.slice().sort((a,b)=>a-b);
  const et = emgTimes.slice().sort((a,b)=>a-b);

  let best = { shift: 0, score: -1, matchPairs: [] };

  // two-pointer match count for a given shift
  function countMatches(shift) {
    let i = 0, j = 0;
    const pairs = [];
    while (i < et.length && j < vt.length) {
      const a = et[i] + shift;
      const b = vt[j];
      const d = a - b;
      if (Math.abs(d) <= tolSec) {
        pairs.push([i, j]);
        i++; j++;
      } else if (d < -tolSec) {
        i++;
      } else {
        j++;
      }
    }
    return pairs;
  }

  for (let s = -rangeSec; s <= rangeSec + 1e-9; s += stepSec) {
    const pairs = countMatches(s);
    const score = pairs.length;
    if (score > best.score) best = { shift: s, score, matchPairs: pairs };
  }

  // recover matched indices in original order (using nearest)
  // We'll use matchPairs' i indices referring to sorted et; to map back we'd need original.
  // Simpler: we will match again on the original list using best.shift.
  return { shift: best.shift, score: best.score };
}

function keepHsMatchedToVideo(hsIdx, tSec, videoEvents, syncRange, tolMs) {
  const tolSec = tolMs / 1000;
  const emgTimes = hsIdx.map((i) => tSec[i]);

  // search shift
  const stepSec = 0.01; // 10 ms
  const { shift, score } = bestShiftMatch(emgTimes, videoEvents, syncRange, stepSec, tolSec);

  // keep matched (two-pointer on sorted)
  const pairs = [];
  const et = emgTimes.map((v, k) => ({ v, k })).sort((a,b)=>a.v-b.v);
  const vt = videoEvents.slice().sort((a,b)=>a-b);

  let i = 0, j = 0;
  while (i < et.length && j < vt.length) {
    const a = et[i].v + shift;
    const b = vt[j];
    const d = a - b;
    if (Math.abs(d) <= tolSec) {
      pairs.push(et[i].k); // index into hsIdx
      i++; j++;
    } else if (d < -tolSec) i++;
    else j++;
  }

  const keepSet = new Set(pairs);
  const hsKept = hsIdx.filter((_, k) => keepSet.has(k));

  return { hsKept, shift, matched: hsKept.length, raw: hsIdx.length, score };
}

/* =======================
   Downloads
   ======================= */

function makeHSCSV(t, hs) {
  let s = "hs_index,No,time_s\n";
  for (let i = 0; i < hs.length; i++) {
    const no = hs[i] + 1;
    const ts = t[hs[i]];
    s += `${i + 1},${no},${ts}\n`;
  }
  return s;
}
function downloadText(filename, text) {
  const blob = new Blob([text], { type: "text/plain;charset=utf-8" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  a.remove();
  URL.revokeObjectURL(url);
}
async function downloadPlotPNG(plotId, filename) {
  const div = document.getElementById(plotId);
  const url = await Plotly.toImage(div, { format: "png", height: 800, width: 1200, scale: 2 });
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  a.remove();
}

/* =======================
   UI & run
   ======================= */

function clearPlots() {
  for (const id of ["plot1", "plot2", "debug1", "debug2"]) {
    const el = document.getElementById(id);
    el.innerHTML = "";
  }
  ["hsn1", "hsn2", "kmad1", "kmad2", "cycn1", "cycn2", "hssrc", "hssrc2"].forEach((id) => ($(id).textContent = "-"));
  $("name1").textContent = "CSV 1（TA推奨）";
  $("name2").textContent = "CSV 2（測定筋：TAのHSで周期化）";
}

async function handleFile(file, sigSel, timeSel, fs) {
  if (!file) return null;

  const { buf, text } = await readFileData(file);
  const parsed = parseCSV(buf, text, file.name, fs);

  if (parsed.format === "table") {
    populateSelect(sigSel, parsed.columns, parsed.guessedSigCol);
    populateSelect(timeSel, parsed.columns, parsed.guessedTimeCol);
  } else {
    populateSelect(sigSel, ["value"], "value");
    populateSelect(timeSel, ["No"], "No");
  }
  return parsed;
}

async function run() {
  try {
    setMsg("読み込み中…");
    clearPlots();

    const fs = parseFloat($("fs").value);
    const hp = parseFloat($("hp").value);
    const lp = parseFloat($("lp").value);
    const rmsMs = parseFloat($("rmsms").value);
    const minBurstMs = parseFloat($("minburst").value);
    const minGapMs = parseFloat($("mingap").value);
    const offsetMs = parseFloat($("offset").value);
    const steps = parseInt($("steps").value, 10);
    const nPoints = parseInt($("npoints").value, 10);
    const showIndividual = $("showind").value === "1";

    const baselineSec = parseFloat($("baselineSec").value);
    const useVideo = $("useVideoAssist").checked && Array.isArray(state.videoEvents) && state.videoEvents.length >= 4;
    const tolMs = parseFloat($("matchTol").value);
    const syncRange = parseFloat($("syncRange").value);

    if (!state.file1 || !state.file2) {
      setMsg("CSVを2つ選んでください。", true);
      return;
    }

    state.parsed1 = await handleFile(state.file1, $("sigcol1"), $("timecol1"), fs);
    state.parsed2 = await handleFile(state.file2, $("sigcol2"), $("timecol2"), fs);

    // expected steps: CSV1(TA) is one side. If total steps known, approx half.
    const expTA = (steps > 0) ? Math.floor(steps / 2) : null;

    setMsg("解析中…");

    const s1 = getSelectedSignal(state.parsed1, $("sigcol1"), $("timecol1"), fs); // TA
    const s2 = getSelectedSignal(state.parsed2, $("sigcol2"), $("timecol2"), fs); // muscle

    const params = { hp, lp, rmsMs, minBurstMs, minGapMs, offsetMs, nPoints };

    // compute independently first (TA HS candidates)
    state.result1 = computePipeline(s1.t, s1.x, fs, params, expTA);
    state.result2 = computePipeline(s2.t, s2.x, fs, params, null); // HS here will be overwritten by TA HS

    // TA HS refinement
    let hsRef = state.result1.hs.slice();

    // remove events before standing baseline
    hsRef = hsRef.filter((i) => state.result1.t[i] >= baselineSec);

    let hsSource = "TA(MAD)";
    if (useVideo) {
      const { hsKept, shift, matched, raw } = keepHsMatchedToVideo(
        hsRef,
        state.result1.t,
        state.videoEvents.filter((v) => v >= baselineSec),
        syncRange,
        tolMs
      );

      if (matched >= 4) {
        hsRef = hsKept;
        hsSource = `TA+Video (shift=${shift.toFixed(2)}s, keep ${matched}/${raw})`;
      } else {
        // fallback if video matching weak
        hsRef = chooseEveryOtherIfDouble(hsRef, fs, expTA);
        hsSource = "TA(1つおき)";
      }
    } else {
      hsRef = chooseEveryOtherIfDouble(hsRef, fs, expTA);
      if (hsRef.length !== state.result1.hs.length) hsSource = "TA(1つおき)";
    }

    // Apply HS reference to TA and muscle
    state.result1 = applyHsToExistingResult(state.result1, hsRef, nPoints);
    state.result2 = applyHsToExistingResult(state.result2, hsRef, nPoints);

    $("name1").textContent = state.parsed1.name || "CSV 1（TA）";
    $("name2").textContent = state.parsed2.name || "CSV 2（測定筋）";

    $("hsn1").textContent = String(state.result1.hs.length);
    $("hsn2").textContent = String(state.result2.hs.length);
    $("kmad1").textContent = String(state.result1.k);
    $("kmad2").textContent = "(TA参照)";
    $("cycn1").textContent = String(Math.max(0, state.result1.hs.length - 1));
    $("cycn2").textContent = String(Math.max(0, state.result2.hs.length - 1));
    $("hssrc").textContent = hsSource;
    $("hssrc2").textContent = hsSource;

    await plotCycle("plot1", state.result1, showIndividual);
    await plotDebug("debug1", state.result1);
    await plotCycle("plot2", state.result2, showIndividual);
    await plotDebug("debug2", state.result2);

    setMsg("完了。TAのHSでCSV2を周期化しています。HS一覧CSV/PNG保存ボタンで出力できます。");
  } catch (e) {
    console.error(e);
    setMsg(String(e?.message || e), true);
  }
}

/* =======================
   Wire UI
   ======================= */

$("file1").addEventListener("change", (ev) => {
  state.file1 = ev.target.files?.[0] || null;
  if (state.file1) setMsg(`CSV1 選択: ${state.file1.name}`);
});
$("file2").addEventListener("change", (ev) => {
  state.file2 = ev.target.files?.[0] || null;
  if (state.file2) setMsg(`CSV2 選択: ${state.file2.name}`);
});

$("run").addEventListener("click", run);

$("clear").addEventListener("click", () => {
  state.file1 = null;
  state.file2 = null;
  $("file1").value = "";
  $("file2").value = "";
  clearPlots();
  setMsg("クリアしました。");
});

// download buttons
$("dlhs1").addEventListener("click", () => {
  if (!state.result1) return;
  const csv = makeHSCSV(state.result1.t, state.result1.hs);
  downloadText(`${state.parsed1?.name || "TA"}_HS.csv`, csv);
});
$("dlhs2").addEventListener("click", () => {
  if (!state.result2) return;
  const csv = makeHSCSV(state.result2.t, state.result2.hs);
  downloadText(`${state.parsed2?.name || "muscle"}_HS_byTA.csv`, csv);
});
$("dlpng1").addEventListener("click", async () => {
  if (!state.result1) return;
  await downloadPlotPNG("plot1", `${state.parsed1?.name || "TA"}_cycle.png`);
});
$("dlpng2").addEventListener("click", async () => {
  if (!state.result2) return;
  await downloadPlotPNG("plot2", `${state.parsed2?.name || "muscle"}_cycle_byTA.png`);
});

// video
$("videoFile").addEventListener("change", (ev) => {
  const file = ev.target.files?.[0] || null;
  state.videoFile = file;
  state.videoEvents = null;
  state.videoEnergy = null;
  state.roi = null;

  const video = $("video");
  const overlay = $("overlay");
  overlay.getContext("2d").clearRect(0,0,overlay.width, overlay.height);

  if (!file) {
    video.removeAttribute("src");
    setVideoStatus("動画未選択");
    return;
  }
  const url = URL.createObjectURL(file);
  video.src = url;
  video.load();
  setVideoStatus(`動画選択: ${file.name}（足付近をタップしてROI指定）`);
});

$("video").addEventListener("loadedmetadata", () => {
  drawRoiOverlay();
  setVideoStatus("動画読込OK：足付近を1回タップしてROI指定 → 「動画からイベント抽出」");
});

// ROI tap (use click/touch coordinates)
$("video").addEventListener("click", (ev) => {
  setRoiFromTap(ev.clientX, ev.clientY);
  setVideoStatus("ROI指定OK（緑枠）。「動画からイベント抽出」を押してください。");
});
window.addEventListener("resize", () => drawRoiOverlay());
$("roiSize").addEventListener("change", () => drawRoiOverlay());

$("analyzeVideo").addEventListener("click", async () => {
  try {
    await extractVideoEvents();
  } catch (e) {
    console.error(e);
    setVideoStatus(String(e?.message || e), true);
  }
});

$("clearVideo").addEventListener("click", () => {
  state.videoFile = null;
  state.videoEvents = null;
  state.videoEnergy = null;
  state.roi = null;

  const video = $("video");
  try {
    if (video.src) URL.revokeObjectURL(video.src);
  } catch {}
  video.removeAttribute("src");
  video.load();
  $("videoFile").value = "";
  $("videoPlot").innerHTML = "";
  drawRoiOverlay();
  setVideoStatus("動画クリアしました。");
});

clearPlots();
setMsg("CSVを2つ選択して「解析してプロット」を押してください。");
setVideoStatus("動画未選択（補助なしでも解析できます）");
