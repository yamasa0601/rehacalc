/**
 * EMG walking app (serverless, browser-only).
 * Pipeline:
 *  - read CSV (sensor-style or table)
 *  - bandpass (HP then LP; 2nd-order Butterworth biquad cascaded twice; forward-backward for ~zero-phase)
 *  - rectify
 *  - RMS envelope (centered)
 *  - HS detection from burst onsets (median + k*MAD), choose k by scoring
 *  - %peak normalization (max=100)
 *  - cycle normalization (0-100%, N points)
 *  - plot mean±SD (+ individual) with Plotly
 */

const $ = (id) => document.getElementById(id);

const state = {
  file1: null,
  file2: null,
  parsed1: null,
  parsed2: null,
  result1: null,
  result2: null,
};

function setMsg(s, isErr = false) {
  const el = $("msg");
  el.textContent = s;
  el.style.color = isErr ? "#ffb4b4" : "var(--muted)";
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

function linspace(a, b, n) {
  const out = new Float64Array(n);
  const step = (b - a) / (n - 1);
  for (let i = 0; i < n; i++) out[i] = a + step * i;
  return out;
}

/** Biquad coefficients using RBJ cookbook. Q=1/sqrt(2) for Butterworth 2nd order */
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
  let x1 = 0,
    x2 = 0,
    y1 = 0,
    y2 = 0;
  for (let i = 0; i < x.length; i++) {
    const x0 = x[i];
    const y0 = c.b0 * x0 + c.b1 * x1 + c.b2 * x2 - c.a1 * y1 - c.a2 * y2;
    y[i] = y0;
    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;
  }
  return y;
}

function reverseArray(x) {
  const y = new Float64Array(x.length);
  for (let i = 0; i < x.length; i++) y[i] = x[x.length - 1 - i];
  return y;
}

/** forward-backward filtering (like filtfilt, approximate) */
function filtfilt(x, coeffs) {
  let y = biquadFilter(x, coeffs);
  y = reverseArray(y);
  y = biquadFilter(y, coeffs);
  y = reverseArray(y);
  return y;
}

/** 4th-order-ish by applying the same 2nd-order twice */
function butterworth2xFILT(x, type, fc, fs) {
  const c = biquadCoeffs(type, fc, fs, 1 / Math.SQRT2);
  let y = filtfilt(x, c);
  y = filtfilt(y, c);
  return y;
}

function demean(x) {
  let s = 0,
    n = 0;
  for (const v of x) {
    if (Number.isFinite(v)) {
      s += v;
      n++;
    }
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

/** centered moving RMS with variable window at edges */
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

function detectHS(env, fs, opts) {
  const {
    expectedCount = null,
    minBurstMs = 30,
    minGapMs = 300,
    offsetMs = 0,
    plausibleStep = [0.35, 1.6],
  } = opts;

  const med = median(env);
  const m = mad(env);
  const minBurst = Math.round((fs * minBurstMs) / 1000);
  const minGap = Math.round((fs * minGapMs) / 1000);
  const offset = Math.round((fs * offsetMs) / 1000);

  const ks = [];
  for (let k = 1.2; k <= 4.0 + 1e-9; k += 0.1) {
    ks.push(Math.round(k * 100) / 100);
  }

  let best = null;

  for (const k of ks) {
    const thr = med + k * m;

    const mask = new Array(env.length);
    for (let i = 0; i < env.length; i++) mask[i] = env[i] > thr;

    let segs = risingEdges(mask);
    segs = segs.filter(([s, e]) => e - s >= minBurst);
    if (segs.length === 0) continue;

    let hs = segs.map(([s, _e]) => s);
    hs.sort((a, b) => a - b);

    // merge close detections
    const merged = [hs[0]];
    for (let i = 1; i < hs.length; i++) {
      if (hs[i] - merged[merged.length - 1] >= minGap) merged.push(hs[i]);
    }
    hs = merged
      .map((v) => v + offset)
      .filter((v) => v >= 0 && v < env.length);

    if (hs.length < 2) continue;

    // score
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

    let score = 5.0 * bad + 2.0 * cv + 0.5 * Math.abs(medInt - 0.7);
    if (expectedCount !== null) score += 0.8 * Math.abs(hs.length - expectedCount);

    if (best === null || score < best.score) {
      best = {
        k,
        thr,
        hs,
        score,
        stats: { n_hs: hs.length, med_step_s: medInt, cv_step: cv, bad_ratio: bad },
      };
    }
  }

  if (best === null) {
    throw new Error(
      "HS推定に失敗しました（信号が弱い/ノイズ多い/閾値設定不適）。minGapや歩数入力を調整してください。"
    );
  }
  return best;
}

function percentPeak(env, hs) {
  let ref = NaN;
  if (hs && hs.length >= 2) {
    const s = hs[0],
      e = hs[hs.length - 1];
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
    const s = hs[i],
      e = hs[i + 1];
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
    let s = 0,
      c = 0;
    for (let i = 0; i < n; i++) {
      const v = cycles[i][j];
      if (Number.isFinite(v)) {
        s += v;
        c++;
      }
    }
    mean[j] = c > 0 ? s / c : NaN;
  }

  for (let j = 0; j < m; j++) {
    let s = 0,
      c = 0;
    for (let i = 0; i < n; i++) {
      const v = cycles[i][j];
      if (Number.isFinite(v)) {
        s += (v - mean[j]) * (v - mean[j]);
        c++;
      }
    }
    sd[j] = c > 1 ? Math.sqrt(s / c) : NaN;
  }

  return { mean, sd };
}

/* =======================
   iPhone Safari対策：CSV読み取り
   ======================= */

// UTF-8が読めるなら読む。読めなくても buf は必ず使えるようにする
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

// 1byte=1char（latin1）として安全に文字列化（TextDecoder iso-8859-1 → fallback）
function bytesToLatin1String(buf) {
  try {
    return new TextDecoder("iso-8859-1", { fatal: false }).decode(buf);
  } catch (e) {
    const bytes = new Uint8Array(buf);
    let out = "";
    const CHUNK = 4096; // Safari対策で小さめ
    for (let i = 0; i < bytes.length; i += CHUNK) {
      const sub = bytes.subarray(i, i + CHUNK);
      out += String.fromCharCode.apply(null, sub);
    }
    return out;
  }
}

// buf（latin1）を主に使ってセンサ形式を確実に拾う。表形式はUTF-8 textがダメならlatin1にfallback
function parseCSV(buf, text, filename, fsDefault) {
  const name = filename.replace(/\.csv$/i, "");

  const raw = bytesToLatin1String(buf);
  const rawLines = raw.split(/\r\n|\n|\r/);

  // センサ形式：どこかの行に "No," がある（前にメタ情報があってOK）
  let start = -1;
  let headerLine = "";
  for (let i = 0; i < Math.min(rawLines.length, 5000); i++) {
    const line = rawLines[i];
    if (!line) continue;
    const idx = line.indexOf("No,");
    if (idx >= 0) {
      start = i;
      headerLine = line.slice(idx); // BOMっぽいゴミがあっても "No," 以降を使う
      break;
    }
  }

  if (start >= 0) {
    const header = headerLine.split(",");
    const idxNo = header.findIndex((h) => (h || "").trim() === "No");

    let idxV = header.findIndex((h) => (h || "").includes("mV") || (h || "").toLowerCase().includes("voltage") || (h || "").includes("電圧"));
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

  // 表形式CSV（UTF-8が読めないならlatin1の先頭行を使う）
  let lines = (text || "").split(/\r\n|\n|\r/);
  if (!lines[0] || lines[0].split(",").length < 2) {
    // UTF-8がダメそう → rawLinesの最初の非空行をheaderに
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
    for (let i = 0; i < rows.length; i++) {
      cols[key][i] = parseFloat(rows[i][j]);
    }
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

  const sig = sigSel.value || parsed.guessedSigCol;
  const tim = timeSel.value || parsed.guessedTimeCol;

  const x = parsed.table[sig];
  const t0 = parsed.table[tim];

  const maxT = Math.max(...t0);
  const estDur = t0.length / fs;

  // time列が秒っぽいならそのまま、インデックスっぽいなら秒に直す
  const isSeconds = Number.isFinite(maxT) && maxT > estDur * 0.8;
  if (isSeconds) {
    return { t: t0, x };
  }

  const tSec = new Float64Array(t0.length);
  for (let i = 0; i < t0.length; i++) tSec[i] = (t0[i] - t0[0]) / fs;
  return { t: tSec, x };
}

function computePipeline(t, x, fs, params, expectedCount) {
  const { hp, lp, rmsMs, minBurstMs, minGapMs, offsetMs, nPoints } = params;

  let y = demean(x);

  if (hp > 0) {
    y = butterworth2xFILT(y, "highpass", hp, fs);
  }
  if (lp > 0 && lp < fs / 2) {
    y = butterworth2xFILT(y, "lowpass", lp, fs);
  }

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
    envPct,
    envRaw: env,
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
    margin: { l: 50, r: 20, t: 30, b: 45 },
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

  const shapes = res.hs.map((idx) => ({
    type: "line",
    x0: res.t[idx],
    x1: res.t[idx],
    y0: 0,
    y1: 1,
    xref: "x",
    yref: "paper",
    line: { width: 1, color: "rgba(255,255,255,0.25)" },
  }));

  const traces = [
    {
      x: t,
      y: y,
      mode: "lines",
      line: { width: 2, color: "#22c55e" },
      name: "%peak env",
    },
  ];

  const layout = {
    margin: { l: 50, r: 20, t: 30, b: 45 },
    xaxis: { title: "time (s)" },
    yaxis: { title: "%RMS (%peak)" },
    shapes,
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
  for (const id of ["plot1", "plot2", "debug1", "debug2"]) {
    const el = document.getElementById(id);
    el.innerHTML = "";
  }
  ["hsn1", "hsn2", "kmad1", "kmad2", "cycn1", "cycn2"].forEach((id) => ($(id).textContent = "-"));
  $("name1").textContent = "CSV 1";
  $("name2").textContent = "CSV 2";
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

    if (!state.file1 || !state.file2) {
      setMsg("CSVを2つ選んでください。", true);
      return;
    }

    state.parsed1 = await handleFile(state.file1, $("sigcol1"), $("timecol1"), fs);
    state.parsed2 = await handleFile(state.file2, $("sigcol2"), $("timecol2"), fs);

    const exp1 = steps > 0 ? Math.floor(steps / 2) : null;
    const exp2 = steps > 0 ? steps - Math.floor(steps / 2) : null;

    setMsg("解析中…（スマホだと数秒〜十数秒かかることがあります）");

    const s1 = getSelectedSignal(state.parsed1, $("sigcol1"), $("timecol1"), fs);
    const s2 = getSelectedSignal(state.parsed2, $("sigcol2"), $("timecol2"), fs);

    const params = { hp, lp, rmsMs, minBurstMs, minGapMs, offsetMs, nPoints };

    state.result1 = computePipeline(s1.t, s1.x, fs, params, exp1);
    state.result2 = computePipeline(s2.t, s2.x, fs, params, exp2);

    $("name1").textContent = state.parsed1.name;
    $("name2").textContent = state.parsed2.name;

    $("hsn1").textContent = String(state.result1.hs.length);
    $("hsn2").textContent = String(state.result2.hs.length);
    $("kmad1").textContent = String(state.result1.k.toFixed(2));
    $("kmad2").textContent = String(state.result2.k.toFixed(2));
    $("cycn1").textContent = String(Math.max(0, state.result1.hs.length - 1));
    $("cycn2").textContent = String(Math.max(0, state.result2.hs.length - 1));

    await plotCycle("plot1", state.result1, showIndividual);
    await plotDebug("debug1", state.result1);
    await plotCycle("plot2", state.result2, showIndividual);
    await plotDebug("debug2", state.result2);

    setMsg("完了。HS一覧CSV/PNG保存ボタンで出力できます。");
  } catch (e) {
    console.error(e);
    setMsg(String(e?.message || e), true);
  }
}

$("file1").addEventListener("change", (ev) => {
  state.file1 = ev.target.files?.[0] || null;
  if (state.file1) setMsg(`CSV 1 選択: ${state.file1.name}`);
});
$("file2").addEventListener("change", (ev) => {
  state.file2 = ev.target.files?.[0] || null;
  if (state.file2) setMsg(`CSV 2 選択: ${state.file2.name}`);
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

$("dlhs1").addEventListener("click", () => {
  if (!state.result1) return;
  const csv = makeHSCSV(state.result1.t, state.result1.hs);
  downloadText(`${state.parsed1?.name || "csv1"}_HS.csv`, csv);
});
$("dlhs2").addEventListener("click", () => {
  if (!state.result2) return;
  const csv = makeHSCSV(state.result2.t, state.result2.hs);
  downloadText(`${state.parsed2?.name || "csv2"}_HS.csv`, csv);
});

$("dlpng1").addEventListener("click", async () => {
  if (!state.result1) return;
  await downloadPlotPNG("plot1", `${state.parsed1?.name || "csv1"}_cycle.png`);
});
$("dlpng2").addEventListener("click", async () => {
  if (!state.result2) return;
  await downloadPlotPNG("plot2", `${state.parsed2?.name || "csv2"}_cycle.png`);
});

clearPlots();
setMsg("CSVを2つ選択して「解析してプロット」を押してください。");

