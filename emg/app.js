/* app.js (module) - EMG gait analysis in-browser
   - Load two CSVs
   - Bandpass (HP-LP), rectification, RMS envelope
   - HS estimation from CSV1 (TA) via burst segmentation -> Peak-based HS (near peak)
   - Use TA HS to define 0–100% gait cycle for BOTH channels
   - %peak normalization (max RMS = 100), mean ± SD
   - Plotly plots + export HS csv / PNG
*/

const $ = (id) => document.getElementById(id);

const state = {
  file1: null, // TA
  file2: null, // target muscle
  parsed1: null,
  parsed2: null,
  res1: null,
  res2: null,
};

// ------------------------ UI helpers ------------------------
function setMsg(s, isErr = false) {
  const el = $("msg");
  if (!el) return;
  el.textContent = s;
  el.style.color = isErr ? "#ffb4b4" : "var(--muted)";
}

// ------------------------ basic stats ------------------------
function median(arr) {
  const a = Array.from(arr).filter(Number.isFinite);
  if (a.length === 0) return NaN;
  a.sort((x, y) => x - y);
  const mid = Math.floor(a.length / 2);
  return a.length % 2 === 0 ? (a[mid - 1] + a[mid]) / 2 : a[mid];
}

function madScalar(arr) {
  const med = median(arr);
  if (!Number.isFinite(med)) return NaN;
  const dev = Array.from(arr, (v) => Math.abs(v - med));
  return median(dev);
}

function linspace(a, b, n) {
  const out = new Float64Array(n);
  const step = (b - a) / (n - 1);
  for (let i = 0; i < n; i++) out[i] = a + step * i;
  return out;
}

// ------------------------ filters (biquad + filtfilt) ------------------------
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

function reverseArray(x) {
  const y = new Float64Array(x.length);
  for (let i = 0; i < x.length; i++) y[i] = x[x.length - 1 - i];
  return y;
}

function filtfilt2(x, coeffs) {
  let y = biquadFilter(x, coeffs);
  y = reverseArray(y);
  y = biquadFilter(y, coeffs);
  y = reverseArray(y);
  return y;
}

// apply twice -> 4th-order-ish (good for EMG)
function butter2x(x, type, fc, fs) {
  const c = biquadCoeffs(type, fc, fs, 1 / Math.SQRT2);
  let y = filtfilt2(x, c);
  y = filtfilt2(y, c);
  return y;
}

// ------------------------ signal ops ------------------------
function demean(x) {
  let s = 0, n = 0;
  for (const v of x) {
    if (Number.isFinite(v)) { s += v; n++; }
  }
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
  for (let i = 0; i < Math.min(starts.length, ends.length); i++) {
    segs.push([starts[i], ends[i]]);
  }
  return segs;
}

// ------------------------ robust CSV reading (iPhone friendly) ------------------------
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

// supports:
//  - sensor-style: includes line with "No," and voltage column
//  - table-style: header row (time, RF/BF/TA/MG etc.)
function parseCSV(buf, text, filename, fsDefault) {
  const name = filename.replace(/\.csv$/i, "");
  const raw = bytesToLatin1String(buf);
  const rawLines = raw.split(/\r\n|\n|\r/);

  // sensor-style: find "No,"
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

  // table-style
  let lines = (text || "").split(/\r\n|\n|\r/);
  if (!lines[0] || lines[0].split(",").length < 2) {
    lines = rawLines.filter((l) => l && l.includes(","));
  }
  const header = (lines[0] || "")
    .split(",")
    .map((s) => (s || "").trim())
    .filter(Boolean);

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

  let timeCol =
    header.find((h) => ["time", "time_s", "t", "sec", "seconds"].includes((h || "").toLowerCase())) || null;
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
  if (!sel) return;
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
  if (parsed.format === "sensor") {
    return { t: parsed.t, x: parsed.x };
  }

  const sig = (sigSel && sigSel.value) || parsed.guessedSigCol;
  const tim = (timeSel && timeSel.value) || parsed.guessedTimeCol;

  const x = parsed.table[sig];
  const t0 = parsed.table[tim];

  // time column seconds-like?
  const maxT = Math.max(...t0);
  const estDur = t0.length / fs;
  const isSeconds = Number.isFinite(maxT) && maxT > estDur * 0.8;

  if (isSeconds) return { t: t0, x };

  const tSec = new Float64Array(t0.length);
  for (let i = 0; i < t0.length; i++) tSec[i] = (t0[i] - t0[0]) / fs;
  return { t: tSec, x };
}

// ------------------------ Peak-based HS estimation ------------------------
function detectHS(env, fs, opts) {
  const {
    expectedCount = null,  // optional (e.g., steps count/2)
    minBurstMs = 30,
    minGapMs = 650,
    offsetMs = 80,

    plausibleStep = [0.6, 2.5], // seconds
    enableEveryOtherFix = true, // auto "skip every other" if double events
    peakBackMs = 180,
    peakForwardMs = 30,
  } = opts;

  const minBurst = Math.max(1, Math.round((minBurstMs / 1000) * fs));
  const minGap = Math.max(1, Math.round((minGapMs / 1000) * fs));
  const offset = Math.round((offsetMs / 1000) * fs);
  const back = Math.round((peakBackMs / 1000) * fs);
  const fwd  = Math.round((peakForwardMs / 1000) * fs);

  const med = median(env);
  let md = madScalar(env);

  // MADが極小なら救済（平坦や読み取り失敗に耐性）
  if (!Number.isFinite(md) || md < 1e-12) {
    const a = Array.from(env).filter(Number.isFinite).sort((x, y) => x - y);
    if (a.length >= 10) {
      const p25 = a[Math.floor(a.length * 0.25)];
      const p75 = a[Math.floor(a.length * 0.75)];
      md = (p75 - p25) / 1.349;
    } else {
      md = 1e-6;
    }
  }

  function scoreHS(hs) {
    if (!hs || hs.length < 2) return Infinity;

    const intervals = [];
    for (let i = 1; i < hs.length; i++) intervals.push((hs[i] - hs[i - 1]) / fs);

    const meanInt = intervals.reduce((a, b) => a + b, 0) / intervals.length;
    const medInt = median(intervals);
    const sdInt = Math.sqrt(
      intervals.reduce((a, b) => a + (b - meanInt) * (b - meanInt), 0) / intervals.length
    );
    const cv = sdInt / (meanInt + 1e-12);

    const [lo, hi] = plausibleStep;
    const bad = intervals.filter((v) => v < lo || v > hi).length / intervals.length;

    let score = 5.0 * bad + 2.0 * cv + 0.5 * Math.abs(medInt - 1.0);
    if (expectedCount !== null) score += 0.8 * Math.abs(hs.length - expectedCount);
    return score;
  }

  function mergeByMinGap(hs, strength) {
    if (hs.length === 0) return { hs, strength };
    const outH = [hs[0]];
    const outS = [strength[0]];
    for (let i = 1; i < hs.length; i++) {
      const last = outH[outH.length - 1];
      if (hs[i] - last < minGap) {
        if (strength[i] > outS[outS.length - 1]) {
          outH[outH.length - 1] = hs[i];
          outS[outS.length - 1] = strength[i];
        }
      } else {
        outH.push(hs[i]);
        outS.push(strength[i]);
      }
    }
    return { hs: outH, strength: outS };
  }

  function chooseEveryOtherIfDouble(hs) {
    if (!enableEveryOtherFix || !hs || hs.length < 6) return hs;

    const intervals = [];
    for (let i = 1; i < hs.length; i++) intervals.push((hs[i] - hs[i - 1]) / fs);
    const medInt = median(intervals);

    const looksDouble =
      (Number.isFinite(medInt) && medInt < 0.50) ||
      (expectedCount !== null && hs.length >= expectedCount * 1.6);

    if (!looksDouble) return hs;

    const candA = hs.filter((_, i) => i % 2 === 0);
    const candB = hs.filter((_, i) => i % 2 === 1);
    return scoreHS(candA) <= scoreHS(candB) ? candA : candB;
  }

  // Peak周辺でHS点を決める：最大傾き点（推奨）
  function hsAroundPeak(segStart, segEnd, peakIdx) {
    const w0 = Math.max(segStart + 1, peakIdx - back);
    const w1 = Math.min(segEnd - 1, peakIdx + fwd);
    let bestI = w0;
    let bestSlope = -Infinity;
    for (let i = w0; i <= w1; i++) {
      const slope = env[i] - env[i - 1];
      if (slope > bestSlope) {
        bestSlope = slope;
        bestI = i;
      }
    }
    return bestI;
  }

  // k sweep
  const ks = [];
  for (let k = 2.0; k <= 6.0 + 1e-12; k += 0.25) ks.push(Number(k.toFixed(2)));

  let best = null;

  for (const k of ks) {
    const thr = med + k * (md * 1.4826); // MAD->sigma
    const mask = new Array(env.length);
    for (let i = 0; i < env.length; i++) mask[i] = env[i] > thr;

    let segs = risingEdges(mask);
    segs = segs
      .map(([s, e]) => [s, e - 1]) // risingEdges gives end as exclusive in our earlier version; normalize
      .map(([s, e]) => [s, Math.max(s, e)])
      .filter(([s, e]) => (e - s + 1) >= minBurst);

    if (segs.length === 0) continue;

    const hs = [];
    const strength = [];

    for (const [s, e] of segs) {
      // find Peak (max amp) in burst
      let p = s;
      let pv = -Infinity;
      for (let i = s; i <= e; i++) {
        const v = env[i];
        if (v > pv) { pv = v; p = i; }
      }
      // HS near peak (max slope) + offset
      const h = hsAroundPeak(s, e, p) + offset;
      if (h >= 0 && h < env.length) {
        hs.push(h);
        strength.push(pv);
      }
    }

    // sort
    const ord = hs.map((v, i) => [v, i]).sort((a, b) => a[0] - b[0]);
    let hsS = ord.map(([v]) => v);
    let stS = ord.map(([, i]) => strength[i]);

    // merge and every-other fix
    const merged = mergeByMinGap(hsS, stS);
    let hsFinal = merged.hs;
    hsFinal = chooseEveryOtherIfDouble(hsFinal);

    if (hsFinal.length < 2) continue;

    const sc = scoreHS(hsFinal);
    if (best === null || sc < best.score) {
      best = { k, thr, hs: hsFinal, score: sc, stats: { n_hs: hsFinal.length } };
    }
  }

  if (!best || best.hs.length < 2) {
    throw new Error("HS推定に失敗しました（周期が作れません）。minGapや歩数入力、閾値を調整してください。");
  }

  // HSが少ない場合は警告（throwしない）
  if (best.hs.length < 3) {
    console.warn("HSが少ないため、平均/SDの信頼性が低いです（推定HS数 < 3）。");
  }

  return best;
}

// ------------------------ %peak + cycle normalization ------------------------
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
  if (!Number.isFinite(ref) || ref <= 0) {
    const out = new Float64Array(env.length);
    out.fill(NaN);
    return out;
  }
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

// ------------------------ pipeline ------------------------
function computePipeline(t, x, fs, params, expectedCount) {
  const {
    hp, lp, rmsMs,
    minBurstMs, minGapMs, offsetMs,
    plausibleStep, enableEveryOtherFix,
    peakBackMs, peakForwardMs,
    nPoints,
  } = params;

  let y = demean(x);
  if (hp > 0) y = butter2x(y, "highpass", hp, fs);
  if (lp > 0 && lp < fs / 2) y = butter2x(y, "lowpass", lp, fs);

  y = absArray(y);
  const win = Math.max(1, Math.round((fs * rmsMs) / 1000));
  const env = movingRMS(y, win);

  const hsRes = detectHS(env, fs, {
    expectedCount,
    minBurstMs,
    minGapMs,
    offsetMs,
    plausibleStep,
    enableEveryOtherFix,
    peakBackMs,
    peakForwardMs,
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
    stats: hsRes.stats,
    grid: cyc.grid,
    cycles: cyc.cycles,
    mean: ms.mean,
    sd: ms.sd,
  };
}

// TA HSで他チャンネルを切り直す
function applyHsToExistingResult(res, hs, nPoints) {
  const hs2 = hs.filter((i) => i >= 0 && i < res.envRaw.length);
  const envPct = percentPeak(res.envRaw, hs2);
  const cyc = normalizeCycles(envPct, hs2, nPoints);
  const ms = meanAndSD(cyc.cycles);
  return { ...res, hs: hs2, envPct, grid: cyc.grid, cycles: cyc.cycles, mean: ms.mean, sd: ms.sd, k: "TA-ref" };
}

// ------------------------ plotting ------------------------
async function plotCycle(targetId, res, showIndividual) {
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
      x, y: lower, mode: "lines", fill: "tonexty",
      fillcolor: "rgba(59,130,246,0.18)", line: { width: 0 },
      hoverinfo: "skip", showlegend: false,
    });
    traces.push({
      x, y: mean, mode: "lines",
      line: { width: 3, color: "#3b82f6" },
      name: "mean",
    });
  }

  const layout = {
    margin: { l: 50, r: 20, t: 30, b: 45 },
    xaxis: { title: "gait cycle (%)", range: [0, 100] },
    yaxis: { title: "%RMS (%peak)" },
    paper_bgcolor: "rgba(0,0,0,0)",
    plot_bgcolor: "rgba(0,0,0,0)",
    font: { color: "#e9eeff" },
  };

  return Plotly.newPlot(targetId, traces, layout, { displayModeBar: false, responsive: true });
}

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

// ------------------------ main run ------------------------
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

function clearPlots() {
  for (const id of ["plot1", "plot2"]) {
    const el = document.getElementById(id);
    if (el) el.innerHTML = "";
  }
  const ids = ["hsn1", "hsn2", "kmad1", "kmad2", "cycn1", "cycn2", "name1", "name2"];
  for (const id of ids) if ($(id)) $(id).textContent = "-";
}

async function run() {
  try {
    setMsg("読み込み中…");
    clearPlots();

    // ---- parameters from UI (fallback defaults if missing) ----
    const fs = parseFloat($("fs")?.value ?? "1000");

    // bandpass (HP-LP). If your UI says "50-450", set hp=50, lp=450.
    const hp = parseFloat($("hp")?.value ?? "50");
    const lp = parseFloat($("lp")?.value ?? "450");

    const rmsMs = parseFloat($("rmsms")?.value ?? "50");
    const minBurstMs = parseFloat($("minburst")?.value ?? "30");
    const minGapMs = parseFloat($("mingap")?.value ?? "650");
    const offsetMs = parseFloat($("offset")?.value ?? "80");
    const nPoints = parseInt($("npoints")?.value ?? "501", 10);
    const showIndividual = ($("showind")?.value ?? "1") === "1";

    // Optional: total steps (video/manual). Used only as expectedCount for HS selection.
    const steps = parseInt($("steps")?.value ?? "0", 10);
    const expectedTA = steps > 0 ? Math.floor(steps / 2) : null;

    if (!state.file1 || !state.file2) {
      setMsg("CSVを2つ選んでください（CSV1=TA, CSV2=測定筋）。", true);
      return;
    }

    state.parsed1 = await handleFile(state.file1, $("sigcol1"), $("timecol1"), fs);
    state.parsed2 = await handleFile(state.file2, $("sigcol2"), $("timecol2"), fs);

    const s1 = getSelectedSignal(state.parsed1, $("sigcol1"), $("timecol1"), fs);
    const s2 = getSelectedSignal(state.parsed2, $("sigcol2"), $("timecol2"), fs);

    setMsg("解析中…");

    const params = {
      hp, lp, rmsMs, minBurstMs, minGapMs, offsetMs,
      plausibleStep: [0.6, 2.5],
      enableEveryOtherFix: true,
      peakBackMs: 180,
      peakForwardMs: 30,
      nPoints,
    };

    // 1) TAでHS推定（Peak方式）
    state.res1 = computePipeline(s1.t, s1.x, fs, params, expectedTA);

    // 2) 測定筋も解析（env等は作る）
    state.res2 = computePipeline(s2.t, s2.x, fs, params, null);

    // 3) 重要：TAのHSで両方の周期(0-100%)を切る
    const hsRef = state.res1.hs;
    state.res1 = applyHsToExistingResult(state.res1, hsRef, nPoints);
    state.res2 = applyHsToExistingResult(state.res2, hsRef, nPoints);

    // UI text
    if ($("name1")) $("name1").textContent = state.parsed1.name;
    if ($("name2")) $("name2").textContent = state.parsed2.name;

    if ($("hsn1")) $("hsn1").textContent = String(state.res1.hs.length);
    if ($("hsn2")) $("hsn2").textContent = String(state.res2.hs.length);
    if ($("kmad1")) $("kmad1").textContent = String(state.res1.k);
    if ($("kmad2")) $("kmad2").textContent = String(state.res2.k);
    if ($("cycn1")) $("cycn1").textContent = String(Math.max(0, state.res1.hs.length - 1));
    if ($("cycn2")) $("cycn2").textContent = String(Math.max(0, state.res2.hs.length - 1));

    await plotCycle("plot1", state.res1, showIndividual);
    await plotCycle("plot2", state.res2, showIndividual);

    setMsg("完了。");
  } catch (e) {
    console.error(e);
    setMsg(String(e?.message || e), true);
  }
}

// ------------------------ wire up ------------------------
$("file1")?.addEventListener("change", (ev) => {
  state.file1 = ev.target.files?.[0] || null;
  if (state.file1) setMsg(`CSV1(TA) 選択: ${state.file1.name}`);
});
$("file2")?.addEventListener("change", (ev) => {
  state.file2 = ev.target.files?.[0] || null;
  if (state.file2) setMsg(`CSV2(測定筋) 選択: ${state.file2.name}`);
});

$("run")?.addEventListener("click", run);
$("clear")?.addEventListener("click", () => {
  state.file1 = null; state.file2 = null;
  if ($("file1")) $("file1").value = "";
  if ($("file2")) $("file2").value = "";
  clearPlots();
  setMsg("クリアしました。");
});

$("dlhs1")?.addEventListener("click", () => {
  if (!state.res1) return;
  const csv = makeHSCSV(state.res1.t, state.res1.hs);
  downloadText(`${state.parsed1?.name || "TA"}_HS.csv`, csv);
});

$("dlhs2")?.addEventListener("click", () => {
  if (!state.res2) return;
  const csv = makeHSCSV(state.res2.t, state.res2.hs);
  downloadText(`${state.parsed2?.name || "MUSCLE"}_HS.csv`, csv);
});

$("dlpng1")?.addEventListener("click", async () => {
  if (!state.res1) return;
  await downloadPlotPNG("plot1", `${state.parsed1?.name || "TA"}_cycle.png`);
});
$("dlpng2")?.addEventListener("click", async () => {
  if (!state.res2) return;
  await downloadPlotPNG("plot2", `${state.parsed2?.name || "MUSCLE"}_cycle.png`);
});

setMsg("CSV1=TA, CSV2=測定筋 を選択して解析してください。");
