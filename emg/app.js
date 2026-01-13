/* app.js (module) - EMG gait analysis in-browser
   - Load two CSVs
   - Bandpass (50-450 Hz), rectification, RMS (50 ms)
   - HS estimation from CSV1 (TA) via burst onset -> HS (offset)
   - Normalize cycles to 0-100%, %peak, mean ± SD
   - Plotly plots + export PNG
   - (Update) Video HS count (rough) is used to keep detected HS count from drifting far.
*/

const $ = (id) => document.getElementById(id);

const state = {
  file1: null,
  file2: null,
  csv1: null,
  csv2: null,
  hsFromTA: null,
  result1: null,
  result2: null,
};

// ---------- utils ----------
function setMsg(msg) {
  $("msg").textContent = msg;
}
function setHSReport(txt) {
  const el = $("hsreport");
  if (!el) return;
  el.textContent = txt || "";
}
function clearPlots() {
  $("plot1").innerHTML = "";
  $("plot2").innerHTML = "";
  $("debug1").innerHTML = "";
  $("debug2").innerHTML = "";
  $("hsn1").textContent = "-";
  $("hsn2").textContent = "-";
  $("kmad1").textContent = "-";
  $("kmad2").textContent = "-";
  $("cycn1").textContent = "-";
  $("cycn2").textContent = "-";
  $("name1").textContent = "CSV 1";
  $("name2").textContent = "CSV 2";
  setMsg("");
  setHSReport("");
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

async function exportPlotPng(divId, filename) {
  const div = $(divId);
  if (!div || !div.data) {
    alert("プロットがありません");
    return;
  }
  const dataUrl = await Plotly.toImage(div, { format: "png", width: 1000, height: 700, scale: 2 });
  const a = document.createElement("a");
  a.href = dataUrl;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  a.remove();
}

// ---------- CSV parsing ----------
function guessDelimiter(text) {
  if (text.includes("\t")) return "\t";
  if (text.includes(";")) return ";";
  return ",";
}

function parseCSV(text) {
  // simple CSV parser (no quotes nesting)
  const delim = guessDelimiter(text);
  const lines = text.split(/\r?\n/).filter((l) => l.trim().length > 0);
  const header = lines[0].split(delim).map((s) => s.trim());
  const rows = [];
  for (let i = 1; i < lines.length; i++) {
    const cols = lines[i].split(delim);
    if (cols.length < header.length) continue;
    const obj = {};
    for (let j = 0; j < header.length; j++) obj[header[j]] = cols[j]?.trim();
    rows.push(obj);
  }
  return { header, rows };
}

function inferColumns(header) {
  const timeKeys = header.filter((h) => /time|t\(s\)|sec|秒/i.test(h));
  const sigKeys = header.filter((h) => !/time|t\(s\)|sec|秒/i.test(h));
  return { timeKeys, sigKeys };
}

function autoDetectNoVoltageEvent(header) {
  const lower = header.map((h) => h.toLowerCase());
  const noIdx = lower.findIndex((h) => h === "no" || h.includes("sample") || h.includes("index"));
  const mvIdx = lower.findIndex((h) => h.includes("mv") || h.includes("emg") || h.includes("voltage"));
  const evIdx = lower.findIndex((h) => h.includes("event"));
  if (noIdx !== -1 && mvIdx !== -1) {
    return { timeKey: header[noIdx], sigKey: header[mvIdx], evKey: evIdx !== -1 ? header[evIdx] : null };
  }
  return null;
}

function buildTimeSignal(rows, timeKey, sigKey, fs) {
  const t = [];
  const x = [];
  for (const r of rows) {
    const tv = Number(r[timeKey]);
    const xv = Number(r[sigKey]);
    if (!Number.isFinite(xv)) continue;

    if (Number.isFinite(tv)) {
      t.push(tv);
    } else {
      const nv = Number(r[timeKey]);
      if (Number.isFinite(nv)) t.push(nv / fs);
      else t.push(t.length / fs);
    }
    x.push(xv);
  }
  const t0 = t[0] ?? 0;
  for (let i = 0; i < t.length; i++) t[i] = t[i] - t0;
  return { t, x };
}

async function readFileText(file) {
  const buf = await file.arrayBuffer();
  try {
    const dec = new TextDecoder("shift_jis");
    const txt = dec.decode(buf);
    if (txt.includes(",") || txt.includes("\n")) return txt;
  } catch {}
  const dec2 = new TextDecoder("utf-8");
  return dec2.decode(buf);
}

async function loadCsvIntoSelectors(file, sigSelId, timeSelId) {
  const txt = await readFileText(file);
  const parsed = parseCSV(txt);
  const auto = autoDetectNoVoltageEvent(parsed.header);

  const sigSel = $(sigSelId);
  const timeSel = $(timeSelId);
  sigSel.innerHTML = "";
  timeSel.innerHTML = "";

  for (const h of parsed.header) {
    const o1 = document.createElement("option");
    o1.value = h; o1.textContent = h;
    sigSel.appendChild(o1);
    const o2 = document.createElement("option");
    o2.value = h; o2.textContent = h;
    timeSel.appendChild(o2);
  }

  if (auto) {
    timeSel.value = auto.timeKey;
    sigSel.value = auto.sigKey;
  } else {
    const { timeKeys, sigKeys } = inferColumns(parsed.header);
    if (timeKeys.length) timeSel.value = timeKeys[0];
    if (sigKeys.length) sigSel.value = sigKeys[0];
  }

  return parsed;
}

// ---------- signal processing ----------
function mean(arr) {
  let s = 0;
  for (const v of arr) s += v;
  return s / arr.length;
}
function std(arr) {
  const m = mean(arr);
  let s = 0;
  for (const v of arr) s += (v - m) * (v - m);
  return Math.sqrt(s / (arr.length - 1));
}

// 2nd order biquad design (RBJ cookbook)
function biquadHPF(fc, fs, Q = 0.7071) {
  const w0 = 2 * Math.PI * (fc / fs);
  const cos = Math.cos(w0);
  const sin = Math.sin(w0);
  const alpha = sin / (2 * Q);

  const b0 = (1 + cos) / 2;
  const b1 = -(1 + cos);
  const b2 = (1 + cos) / 2;
  const a0 = 1 + alpha;
  const a1 = -2 * cos;
  const a2 = 1 - alpha;
  return { b0: b0 / a0, b1: b1 / a0, b2: b2 / a0, a1: a1 / a0, a2: a2 / a0 };
}
function biquadLPF(fc, fs, Q = 0.7071) {
  const w0 = 2 * Math.PI * (fc / fs);
  const cos = Math.cos(w0);
  const sin = Math.sin(w0);
  const alpha = sin / (2 * Q);

  const b0 = (1 - cos) / 2;
  const b1 = 1 - cos;
  const b2 = (1 - cos) / 2;
  const a0 = 1 + alpha;
  const a1 = -2 * cos;
  const a2 = 1 - alpha;
  return { b0: b0 / a0, b1: b1 / a0, b2: b2 / a0, a1: a1 / a0, a2: a2 / a0 };
}

function applyBiquad(x, c) {
  const y = new Array(x.length);
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
function filtfiltBiquad(x, c) {
  const fwd = applyBiquad(x, c);
  const rev = applyBiquad(fwd.slice().reverse(), c).reverse();
  return rev;
}

function bandpassRectifyRms(x, fs, hp, lp, rmsMs) {
  let y = x.slice();
  if (hp && hp > 0) y = filtfiltBiquad(y, biquadHPF(hp, fs));
  if (lp && lp > 0) y = filtfiltBiquad(y, biquadLPF(lp, fs));
  y = y.map((v) => Math.abs(v));
  const w = Math.max(1, Math.round((rmsMs / 1000) * fs));
  const out = new Array(y.length);
  let sumsq = 0;
  const q = [];
  for (let i = 0; i < y.length; i++) {
    const v = y[i];
    sumsq += v * v;
    q.push(v);
    if (q.length > w) {
      const old = q.shift();
      sumsq -= old * old;
    }
    out[i] = Math.sqrt(sumsq / q.length);
  }
  return out;
}

function toPercentPeak(env) {
  const peak = Math.max(...env);
  if (!Number.isFinite(peak) || peak <= 0) return env.map(() => 0);
  return env.map((v) => (v / peak) * 100);
}

// ---------- HS estimation helpers ----------
function mad(arr) {
  const m = arr.slice().sort((a, b) => a - b)[Math.floor(arr.length / 2)];
  const dev = arr.map((v) => Math.abs(v - m)).sort((a, b) => a - b);
  const md = dev[Math.floor(dev.length / 2)];
  return { median: m, mad: md };
}
function pickCandidatesAbove(env, thr) {
  const idx = [];
  for (let i = 0; i < env.length; i++) if (env[i] >= thr) idx.push(i);
  return idx;
}
function groupConsecutive(indices) {
  const groups = [];
  if (indices.length === 0) return groups;
  let s = indices[0], prev = indices[0];
  for (let i = 1; i < indices.length; i++) {
    const v = indices[i];
    if (v === prev + 1) prev = v;
    else { groups.push([s, prev]); s = v; prev = v; }
  }
  groups.push([s, prev]);
  return groups;
}
function hsIntervals(hs, fs) {
  const intervals = [];
  for (let i = 1; i < hs.length; i++) intervals.push((hs[i] - hs[i - 1]) / fs);
  return intervals;
}

// ---------- 任意：動画HS（手入力：一覧 or HS数） ----------
function secToIndex(sec, t) {
  if (!Number.isFinite(sec) || !t || t.length === 0) return NaN;
  const t0 = t[0];
  const tN = t[t.length - 1];
  if (!Number.isFinite(t0) || !Number.isFinite(tN)) return NaN;
  if (sec <= t0) return 0;
  if (sec >= tN) return t.length - 1;

  let lo = 0, hi = t.length - 1;
  while (hi - lo > 1) {
    const mid = (lo + hi) >> 1;
    if (t[mid] < sec) lo = mid;
    else hi = mid;
  }
  return Math.abs(t[lo] - sec) <= Math.abs(t[hi] - sec) ? lo : hi;
}

/**
 * 手入力HSの解釈：
 *  - 1つだけ：HS数（expectedCount）として扱う（"10" / "10HS" / "HS10" もOK）
 *  - 複数：HS一覧として扱う（サンプルNo(1始まり) / 秒 / ms / mm:ss）
 */
function parseManualHS(text, fs, t, nSamples) {
  const raw = (text || "").trim();
  if (!raw) return { mode: "none", hs: null, expectedCount: null, note: "" };

  const tokens = raw
    .replace(/[,;]+/g, " ")
    .split(/\s+/)
    .map((s) => s.trim())
    .filter(Boolean);

  // 1つだけなら「HS数」扱い（数値部分を抽出）
  if (tokens.length === 1) {
    const m = String(tokens[0]).match(/(\d+)/);
    if (m) {
      const v = Number(m[1]);
      if (Number.isFinite(v) && Math.floor(v) === v && v > 0 && v <= 500) {
        return { mode: "count", hs: null, expectedCount: v, note: `手入力HS数=${v}（目視情報として使用）` };
      }
    }
  }

  // 0が含まれていれば0始まり、そうでなければ1始まり優先
  let base0 = false;
  for (const tok of tokens) {
    const m = String(tok).match(/^\d+$/);
    if (m) {
      const vv = Number(tok);
      if (vv === 0) base0 = true;
    }
  }

  const hs = [];
  const durS = t && t.length ? (t[t.length - 1] - t[0]) : (nSamples - 1) / fs;

  for (let tok of tokens) {
    let s = String(tok).trim().toLowerCase();
    if (!s) continue;

    // mm:ss(.sss)
    if (s.includes(":")) {
      const parts = s.split(":");
      if (parts.length === 2) {
        const mm = Number(parts[0]);
        const ss = Number(parts[1]);
        if (Number.isFinite(mm) && Number.isFinite(ss)) {
          const sec = mm * 60 + ss;
          const idx = secToIndex(sec + (t ? t[0] : 0), t);
          if (Number.isFinite(idx)) hs.push(idx);
        }
      }
      continue;
    }

    // unit suffix
    let unit = "";
    if (s.endsWith("ms")) { unit = "ms"; s = s.slice(0, -2); }
    else if (s.endsWith("s")) { unit = "s"; s = s.slice(0, -1); }

    const v = Number(s);
    if (!Number.isFinite(v)) continue;

    const isInteger = Math.floor(v) === v;
    if (unit === "" && isInteger && v >= 0 && v <= (nSamples + 1)) {
      const idx = base0 ? v : (v === 0 ? 0 : v - 1);
      if (idx >= 0 && idx < nSamples) hs.push(idx);
      continue;
    }

    let sec = NaN;
    if (unit === "ms") sec = v / 1000;
    else if (unit === "s") sec = v;
    else {
      if (v > durS && v <= durS * 1000 + 2000) sec = v / 1000;
      else sec = v;
    }

    const idx = t ? secToIndex(sec, t) : Math.round(sec * fs);
    if (Number.isFinite(idx) && idx >= 0 && idx < nSamples) hs.push(idx);
  }

  const hsSorted = Array.from(new Set(hs.filter((x) => Number.isFinite(x)))).sort((a, b) => a - b);
  if (hsSorted.length < 2) {
    return { mode: "bad", hs: null, expectedCount: null, note: "手入力HSの解釈に失敗しました（値が少なすぎる/範囲外）。" };
  }
  return { mode: "list", hs: hsSorted, expectedCount: hsSorted.length, note: `手入力HS一覧=${hsSorted.length}点` };
}

function compareHS(hsAuto, hsManual, fs, tolMs) {
  const tol = Math.max(1, Math.round((fs * tolMs) / 1000));
  const a = (hsAuto || []).slice().sort((x, y) => x - y);
  const m = (hsManual || []).slice().sort((x, y) => x - y);
  let i = 0, j = 0;
  const diffs = [];
  while (i < a.length && j < m.length) {
    const d = a[i] - m[j];
    if (Math.abs(d) <= tol) { diffs.push(d); i++; j++; continue; }
    if (a[i] < m[j] - tol) i++;
    else j++;
  }
  const matched = diffs.length;
  const miss = m.length - matched;
  const extra = a.length - matched;

  const absMs = diffs.map((d) => (Math.abs(d) / fs) * 1000);
  const biasMs = diffs.map((d) => (d / fs) * 1000);
  const mae = absMs.length ? absMs.reduce((s, v) => s + v, 0) / absMs.length : NaN;
  const bias = biasMs.length ? biasMs.reduce((s, v) => s + v, 0) / biasMs.length : NaN;

  return { matched, miss, extra, mae, bias };
}

// (Update) general downsample (every m-th) to match expectedCount
function chooseStrideSubseq(hs, fs, expectedCount, plausibleStep, countTol) {
  if (!expectedCount || expectedCount <= 0) return null;
  if (hs.length < expectedCount + Math.max(1, countTol)) return null;

  const ratio = hs.length / expectedCount;
  const mCandidates = [];
  if (ratio >= 1.6) mCandidates.push(2);
  if (ratio >= 2.6) mCandidates.push(3);

  function score(seq) {
    if (seq.length < 3) return Infinity;
    const ints = hsIntervals(seq, fs);
    const med = ints.slice().sort((a, b) => a - b)[Math.floor(ints.length / 2)];
    const sd = std(ints);
    let pen = 0;
    if (plausibleStep) {
      const [lo, hi] = plausibleStep;
      if (med < lo) pen += (lo - med) * 3;
      if (med > hi) pen += (med - hi) * 3;
    }
    const diff = Math.abs(seq.length - expectedCount);
    // HARD-ish penalty if far
    pen += diff > countTol ? 100 + diff * 10 : diff * 0.5;
    return sd + pen;
  }

  let best = null;
  let bestScore = Infinity;

  for (const m of mCandidates) {
    for (let off = 0; off < m; off++) {
      const seq = hs.filter((_, i) => (i % m) === off);
      const sc = score(seq);
      if (sc < bestScore) {
        bestScore = sc;
        best = { name: `every${m}@${off}`, hs: seq };
      }
    }
  }

  // return only if it improves count proximity
  if (!best) return null;
  const beforeDiff = Math.abs(hs.length - expectedCount);
  const afterDiff = Math.abs(best.hs.length - expectedCount);
  if (afterDiff < beforeDiff) return best;
  return null;
}

function chooseOddEven(hs, fs, expectedCount, plausibleStep, countTol) {
  const subseqs = [];
  if (hs.length >= 4) {
    const odd = hs.filter((_, i) => i % 2 === 0);
    const even = hs.filter((_, i) => i % 2 === 1);
    if (odd.length >= 3) subseqs.push({ name: "odd", hs: odd });
    if (even.length >= 3) subseqs.push({ name: "even", hs: even });
  }
  if (subseqs.length === 0) return null;

  function score(seq) {
    const ints = hsIntervals(seq, fs);
    const med = ints.slice().sort((a, b) => a - b)[Math.floor(ints.length / 2)];
    const sd = std(ints);
    let pen = 0;
    if (plausibleStep) {
      const [lo, hi] = plausibleStep;
      if (med < lo) pen += (lo - med) * 2;
      if (med > hi) pen += (med - hi) * 2;
    }
    if (expectedCount && expectedCount > 0) {
      const diff = Math.abs(seq.length - expectedCount);
      pen += diff > countTol ? 60 + diff * 8 : diff * 0.6;
    }
    return sd + pen;
  }

  let best = subseqs[0];
  let bestScore = score(best.hs);
  for (let i = 1; i < subseqs.length; i++) {
    const sc = score(subseqs[i].hs);
    if (sc < bestScore) { best = subseqs[i]; bestScore = sc; }
  }
  return best;
}

function detectHS(env, fs, opts) {
  const {
    expectedCount = null,
    countTol = 2,
    baselineSec = 3,
    minBurstMs = 30,
    minGapMs = 650,
    offsetMs = 80,
    plausibleStep = [0.6, 2.5],
    targetStepS = 1.0,
    minStrideS = 0.65,
  } = opts;

  const baseN = Math.min(env.length, Math.max(10, Math.round(baselineSec * fs)));
  const base = env.slice(0, baseN);
  const { median, mad: md } = mad(base);
  const sigma = md * 1.4826;

  const ks = [];
  for (let k = 2.0; k <= 6.0; k += 0.25) ks.push(Number(k.toFixed(2)));

  const minLen = Math.max(1, Math.round((minBurstMs / 1000) * fs));
  const minGap = Math.max(1, Math.round((minGapMs / 1000) * fs));
  const offset = Math.round((offsetMs / 1000) * fs);

  function buildHSforK(k) {
    const thr = median + k * sigma;
    const idx = pickCandidatesAbove(env, thr);
    const groups = groupConsecutive(idx).filter(([s, e]) => (e - s + 1) >= minLen);

    let hs = groups.map(([s, _e]) => s + offset).filter((i) => i >= 0 && i < env.length);

    // suppress duplicates within minGap
    const filtered = [];
    for (const h of hs) {
      if (filtered.length === 0 || (h - filtered[filtered.length - 1]) >= minGap) filtered.push(h);
    }
    hs = filtered;

    // 1) odd/even
    const choice = chooseOddEven(hs, fs, expectedCount, plausibleStep, countTol);
    if (choice) hs = choice.hs;

    // 2) more aggressive downsample (2x, 3x etc) when count is much larger than expected
    const strideChoice = chooseStrideSubseq(hs, fs, expectedCount, plausibleStep, countTol);
    if (strideChoice) hs = strideChoice.hs;

    return { k, thr, hs };
  }

  function scoreHS(hs) {
    if (hs.length < 3) return Infinity;

    const ints = hsIntervals(hs, fs);
    const med = ints.slice().sort((a, b) => a - b)[Math.floor(ints.length / 2)];
    const sd = std(ints);

    let pen = 0;
    const [lo, hi] = plausibleStep || [0, Infinity];
    if (med < lo) pen += (lo - med) * 3;
    if (med > hi) pen += (med - hi) * 3;
    pen += Math.abs(med - targetStepS) * 0.5;

    // (Update) count constraint: avoid drifting far from the rough count
    if (expectedCount && expectedCount > 0) {
      const diff = Math.abs(hs.length - expectedCount);
      // if outside tolerance band => huge penalty (almost a hard constraint)
      pen += diff > countTol ? 200 + diff * 20 : diff * 1.0;
    }

    if (med < minStrideS) pen += (minStrideS - med) * 5;
    return sd + pen;
  }

  let best = null;
  let bestScore = Infinity;

  for (const k of ks) {
    const { thr, hs } = buildHSforK(k);
    const sc = scoreHS(hs);
    if (sc < bestScore) {
      bestScore = sc;
      best = { k, thr, hs };
    }
  }

  if (!best || best.hs.length < 3) {
    throw new Error("HS推定に失敗しました（信号が弱い/ノイズ多い/閾値設定不適）。minGapやHS数入力を調整してください。");
  }

  const intervals = hsIntervals(best.hs, fs);
  const med = intervals.slice().sort((a, b) => a - b)[Math.floor(intervals.length / 2)];
  const countNote = (expectedCount && expectedCount > 0)
    ? ` / count=${best.hs.length} (expected≈${expectedCount}±${countTol})`
    : ` / count=${best.hs.length}`;
  const note = `k(MAD)=${best.k.toFixed(2)}${countNote} / median周期=${med.toFixed(2)}s`;

  return { ...best, note, stats: { n_hs: best.hs.length, med_step_s: med } };
}

// ---------- cycle normalization ----------
function interpolateToGrid(x, y, xGrid) {
  const out = new Array(xGrid.length);
  let j = 0;
  for (let i = 0; i < xGrid.length; i++) {
    const xi = xGrid[i];
    while (j < x.length - 2 && x[j + 1] < xi) j++;
    const x0 = x[j], x1 = x[j + 1];
    const y0 = y[j], y1 = y[j + 1];
    const t = (xi - x0) / (x1 - x0);
    out[i] = y0 + t * (y1 - y0);
  }
  return out;
}

function normalizeCycles(t, envPct, hsIdx, nPoints = 501) {
  const cycles = [];
  for (let i = 0; i < hsIdx.length - 1; i++) {
    const s = hsIdx[i], e = hsIdx[i + 1];
    if (e <= s + 5) continue;
    const segT = t.slice(s, e);
    const segY = envPct.slice(s, e);
    const dur = segT[segT.length - 1] - segT[0];
    if (dur <= 0) continue;
    const p = segT.map((tv) => ((tv - segT[0]) / dur) * 100);
    const grid = [];
    for (let k = 0; k < nPoints; k++) grid.push((k / (nPoints - 1)) * 100);
    const yGrid = interpolateToGrid(p, segY, grid);
    cycles.push({ grid, y: yGrid });
  }
  return cycles;
}

function meanStdOfCycles(cycles) {
  const n = cycles.length;
  if (n === 0) return null;
  const m = new Array(cycles[0].y.length).fill(0);
  for (const c of cycles) for (let i = 0; i < m.length; i++) m[i] += c.y[i];
  for (let i = 0; i < m.length; i++) m[i] /= n;

  const s = new Array(m.length).fill(0);
  for (const c of cycles) {
    for (let i = 0; i < m.length; i++) {
      const d = c.y[i] - m[i];
      s[i] += d * d;
    }
  }
  for (let i = 0; i < s.length; i++) s[i] = Math.sqrt(s[i] / Math.max(1, n - 1));
  return { mean: m, sd: s };
}

function applyHS(t, env, hsIdx, nPoints) {
  const pct = toPercentPeak(env);
  const cycles = normalizeCycles(t, pct, hsIdx, nPoints);
  const ms = meanStdOfCycles(cycles);
  return { pct, hs: hsIdx, cycles, ms };
}

// ---------- plotting ----------
function plotCycles(divId, res, title, showIndividual = true) {
  const grid = res.cycles.length ? res.cycles[0].grid : [];
  const traces = [];

  if (showIndividual) {
    for (const c of res.cycles) {
      traces.push({
        x: c.grid, y: c.y,
        mode: "lines",
        line: { width: 1 },
        opacity: 0.25,
        hoverinfo: "skip",
        name: "cycle",
        showlegend: false,
      });
    }
  }

  if (res.ms) {
    traces.push({ x: grid, y: res.ms.mean, mode: "lines", line: { width: 4 }, name: "mean" });
    const upper = res.ms.mean.map((v, i) => v + res.ms.sd[i]);
    const lower = res.ms.mean.map((v, i) => v - res.ms.sd[i]);
    traces.push({ x: grid, y: upper, mode: "lines", line: { width: 0 }, showlegend: false, hoverinfo: "skip" });
    traces.push({ x: grid, y: lower, mode: "lines", fill: "tonexty", line: { width: 0 }, opacity: 0.2, showlegend: false, hoverinfo: "skip" });
  }

  const layout = {
    title: { text: title, font: { size: 16 } },
    xaxis: { title: "gait cycle (%)", range: [0, 100], zeroline: false, gridcolor: "rgba(255,255,255,0.15)" },
    yaxis: { title: "%RMS (%peak)", rangemode: "tozero", zeroline: false, gridcolor: "rgba(255,255,255,0.15)" },
    paper_bgcolor: "rgba(0,0,0,0)",
    plot_bgcolor: "rgba(0,0,0,0)",
    font: { color: "#e9eeff" },
    margin: { l: 60, r: 20, t: 50, b: 50 },
    legend: { x: 1, y: 1 },
  };

  Plotly.newPlot($(divId), traces, layout, { responsive: true, displayModeBar: true, modeBarButtonsToRemove: ["select2d", "lasso2d"] });
}

function plotDebug(divId, t, pct, hs) {
  const traces = [{ x: t, y: pct, mode: "lines", name: "%peak envelope", line: { width: 2 } }];
  if (hs && hs.length) {
    traces.push({ x: hs.map((i) => t[i]), y: hs.map((i) => pct[i]), mode: "markers", name: "HS", marker: { size: 8 } });
  }
  const layout = {
    xaxis: { title: "time (s)", gridcolor: "rgba(255,255,255,0.12)" },
    yaxis: { title: "%peak", gridcolor: "rgba(255,255,255,0.12)" },
    paper_bgcolor: "rgba(0,0,0,0)",
    plot_bgcolor: "rgba(0,0,0,0)",
    font: { color: "#e9eeff" },
    margin: { l: 60, r: 20, t: 20, b: 40 },
    showlegend: true,
  };
  Plotly.newPlot($(divId), traces, layout, { responsive: true });
}

// ---------- main ----------
async function getSelectedSignal(parsed, sigKey, timeKey, fs) {
  return buildTimeSignal(parsed.rows, timeKey, sigKey, fs);
}

async function run() {
  try {
    clearPlots();
    setMsg("解析中...");

    const fs = Number($("fs").value);
    const hp = Number($("hp").value);
    const lp = Number($("lp").value);
    const rmsMs = Number($("rmsms").value);
    const baselineSec = Number($("baseline").value);
    const minBurstMs = Number($("minburst").value);
    const minGapMs = Number($("mingap").value);
    const offsetMs = Number($("offset").value);

    const steps = Number($("steps").value);
    const stepsTol = Number($("stepstol").value);

    const nPoints = Number($("npoints").value);
    const showInd = $("showind").value === "1";

    if (!state.csv1 || !state.csv2) {
      setMsg("CSVを2つ選択してください。");
      return;
    }

    const s1 = await getSelectedSignal(state.csv1, $("sigcol1").value, $("timecol1").value, fs);
    const s2 = await getSelectedSignal(state.csv2, $("sigcol2").value, $("timecol2").value, fs);

    $("name1").textContent = state.file1?.name || "CSV 1";
    $("name2").textContent = state.file2?.name || "CSV 2";

    const envTA = bandpassRectifyRms(s1.x, fs, hp, lp, rmsMs);
    const envM = bandpassRectifyRms(s2.x, fs, hp, lp, rmsMs);

    const manual = parseManualHS($("hsmanual").value, fs, s1.t, envTA.length);
    const tolMs = parseFloat($("hstol").value);
    const hsUseMode = $("hsuse").value;

    // expectedCount priority: manual count > steps > none
    const expectedCount = manual.expectedCount !== null ? manual.expectedCount : (steps > 0 ? steps : null);
    const countTol = Number.isFinite(stepsTol) ? Math.max(0, Math.round(stepsTol)) : 2;

    let hsRes = null;
    let hsUsed = null;
    let hsUsedNote = "";

    try {
      hsRes = detectHS(envTA, fs, {
        expectedCount,
        countTol,
        baselineSec,
        minBurstMs,
        minGapMs,
        offsetMs,
        plausibleStep: [0.6, 2.5],
        targetStepS: 1.0,
        minStrideS: 0.65,
      });
      hsUsed = hsRes.hs;
      hsUsedNote = hsRes.note || "";
    } catch (err) {
      // if auto fails but manual list exists => still proceed
      if (manual.mode === "list" && manual.hs && manual.hs.length >= 3) {
        hsRes = {
          k: NaN,
          thr: NaN,
          hs: manual.hs,
          note: "自動HS失敗→手入力HSで周期分割",
          stats: { n_hs: manual.hs.length, med_step_s: NaN },
        };
        hsUsed = manual.hs;
        hsUsedNote = hsRes.note;
      } else {
        throw err;
      }
    }

    // report
    if (manual.mode === "list" && manual.hs && hsRes && hsRes.hs) {
      const cmp = compareHS(hsRes.hs, manual.hs, fs, Number.isFinite(tolMs) ? tolMs : 80);
      const mae = Number.isFinite(cmp.mae) ? cmp.mae.toFixed(1) : "-";
      const bias = Number.isFinite(cmp.bias) ? cmp.bias.toFixed(1) : "-";
      const msg = `動画HS(${manual.hs.length}) vs 自動HS(${hsRes.hs.length}): 一致${cmp.matched}  不一致(欠損${cmp.miss}/余分${cmp.extra})  MAE=${mae}ms  bias=${bias}ms`;
      setHSReport(msg + (hsUseMode === "use" ? " / 周期分割=手入力HS" : ""));
    } else if ((manual.mode === "count" && manual.expectedCount !== null) || (steps > 0)) {
      const exp = expectedCount ?? "-";
      setHSReport(`目視HS数≈${exp}±${countTol} / 自動HS=${hsRes?.hs?.length ?? "-"}`);
    } else if (manual.mode === "bad") {
      setHSReport(manual.note);
    } else {
      setHSReport("");
    }

    // use manual HS list for segmentation if requested
    if (manual.mode === "list" && manual.hs && hsUseMode === "use") {
      hsUsed = manual.hs;
      hsUsedNote = (hsUsedNote ? hsUsedNote + " / " : "") + "手入力HSで周期分割";
    }

    state.hsFromTA = hsRes;

    state.result1 = applyHS(s1.t, envTA, hsUsed, nPoints);
    state.result2 = applyHS(s2.t, envM, hsUsed, nPoints);
    state.result1.k = hsRes.k;
    state.result2.k = hsRes.k;
    state.result1.thr = hsRes.thr;
    state.result2.thr = hsRes.thr;
    state.result1.note = hsUsedNote;
    state.result2.note = hsUsedNote;

    $("hsn1").textContent = String(state.result1.hs.length);
    $("hsn2").textContent = String(state.result2.hs.length);
    $("cycn1").textContent = String(Math.max(0, state.result1.hs.length - 1));
    $("cycn2").textContent = String(Math.max(0, state.result2.hs.length - 1));
    $("kmad1").textContent = Number.isFinite(state.hsFromTA?.k) ? state.hsFromTA.k.toFixed(2) : "-";
    $("kmad2").textContent = Number.isFinite(state.hsFromTA?.k) ? state.hsFromTA.k.toFixed(2) : "-";

    plotCycles("plot1", state.result1, `${$("name1").textContent}`, showInd);
    plotCycles("plot2", state.result2, `${$("name2").textContent}`, showInd);
    plotDebug("debug1", s1.t, state.result1.pct, state.result1.hs);
    plotDebug("debug2", s2.t, state.result2.pct, state.result2.hs);

    const note = state.result1?.note || "";
    setMsg(note ? `完了。\n(${note})\nHS一覧CSV/PNG保存ボタンで出力できます。` : "完了。HS一覧CSV/PNG保存ボタンで出力できます。");
  } catch (e) {
    console.error(e);
    setMsg(String(e?.message || e));
  }
}

$("file1").addEventListener("change", async (ev) => {
  state.file1 = ev.target.files?.[0] || null;
  if (!state.file1) return;
  state.csv1 = await loadCsvIntoSelectors(state.file1, "sigcol1", "timecol1");
  setMsg(`CSV1 読み込み: ${state.file1.name}`);
});

$("file2").addEventListener("change", async (ev) => {
  state.file2 = ev.target.files?.[0] || null;
  if (!state.file2) return;
  state.csv2 = await loadCsvIntoSelectors(state.file2, "sigcol2", "timecol2");
  setMsg(`CSV2 読み込み: ${state.file2.name}`);
});

$("run").addEventListener("click", run);

$("clear").addEventListener("click", () => {
  state.file1 = null;
  state.file2 = null;
  state.csv1 = null;
  state.csv2 = null;
  state.hsFromTA = null;
  state.result1 = null;
  state.result2 = null;

  $("file1").value = "";
  $("file2").value = "";
  $("hsmanual").value = "";
  $("hsuse").value = "check";
  clearPlots();
});

$("dlhs1").addEventListener("click", () => {
  if (!state.result1?.hs?.length) return alert("HSがありません");
  const fs = Number($("fs").value);
  const lines = ["hs_index,hs_time_s"];
  for (const idx of state.result1.hs) lines.push(`${idx},${(idx / fs).toFixed(6)}`);
  downloadText(`${$("name1").textContent}_HS.csv`, lines.join("\n"));
});

$("dlhs2").addEventListener("click", () => {
  if (!state.result2?.hs?.length) return alert("HSがありません");
  const fs = Number($("fs").value);
  const lines = ["hs_index,hs_time_s"];
  for (const idx of state.result2.hs) lines.push(`${idx},${(idx / fs).toFixed(6)}`);
  downloadText(`${$("name2").textContent}_HS.csv`, lines.join("\n"));
});

$("dlpng1").addEventListener("click", () => exportPlotPng("plot1", `${$("name1").textContent}.png`));
$("dlpng2").addEventListener("click", () => exportPlotPng("plot2", `${$("name2").textContent}.png`));
