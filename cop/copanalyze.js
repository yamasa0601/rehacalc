/* COP analysis in-browser (GitHub Pages friendly)
   Implements the same core method as the provided Python notebook:
   - Trim first/last N seconds
   - Centering
   - 4th order Butterworth lowpass at cutoff (default 10Hz) with filtfilt-like padding
   - Unit conversion to cm
   - 95% confidence ellipse area using chi2.ppf(0.95, df=2)
   - Convex hull area
   - Welch PSD (Hann, noverlap=50%, detrend=mean) and band ratios within 0–3 Hz
*/

(() => {
  const $ = (id) => document.getElementById(id);

  // --- UI ---
  const fileEl = $("file");
  const runBtn = $("run");
  const msgEl = $("msg");
  const metricsEl = $("metrics");
  const dlMetricsBtn = $("dlMetrics");
  const dlBandsBtn = $("dlBands");
  const dlPngBtn = $("dlPng");
  const detailsEl = $("details");

  let lastResult = null;
  let lastFilename = "cop";

  fileEl.addEventListener("change", () => {
    runBtn.disabled = !fileEl.files || !fileEl.files.length;
    resetOutputs();
  });

  runBtn.addEventListener("click", async () => {
    resetOutputs();
    const f = fileEl.files?.[0];
    if (!f) return;

    const fs = parseFloat($("fs").value || "100");
    const cutoff = parseFloat($("cutoff").value || "10");
    const trimSecInput = parseFloat($("trimSec").value || "5");
    const autoTrim = $("autoTrim") ? $("autoTrim").checked : true;
    const resample = $("resample") ? $("resample").checked : true;
    const unit = $("unit").value; // mm or cm
    const unitToCm = (unit === "mm") ? 0.1 : 1.0;

    lastFilename = (f.name || "cop").replace(/\.[^.]+$/, "");

    showMsg("読み込み中…", "warn");

    try {
      const rows = await parseCsvFile(f);
      const parsed = extractColumns(rows);

      if (!parsed.ok) {
        showMsg(parsed.error, "warn");
        return;
      }

      let { time, copx, copy } = parsed;
// Sort by time and de-duplicate timestamps
({ time, copx, copy } = sortByTime(time, copx, copy));

// Normalize time units (seconds / ms / sample index)
const normTime = normalizeTimeUnits(time, fs);
time = normTime.timeSec;
let preprocessNotes = [];
if (normTime.note) preprocessNotes.push(normTime.note);

if (time.length < 20) {
  showMsg("データ点が少なすぎます（20点以上を推奨）。", "warn");
  return;
}

// Optional: resample to uniform grid BEFORE trimming (fills sampling dropouts)
if (resample && time.length >= 2) {
  const rs = resampleLinearUniform(time, copx, copy, fs);
  time = rs.time; copx = rs.copx; copy = rs.copy;
  preprocessNotes.push("欠損補完：等間隔リサンプル（線形補間）を実施しました。");
}
// Trim by seconds using time column (auto adjust when trial is short)
      const fullDur = time[time.length - 1] - time[0];
      let trimSec = trimSecInput;
      let note = preprocessNotes.join(" ");
      const minKeepSec = Math.max(2.0, 20 / fs); // keep at least 2s (and 20 points)
      if (autoTrim && Number.isFinite(fullDur) && fullDur > 0) {
        const maxTrim = Math.max(0, (fullDur - minKeepSec) / 2);
        if (trimSec > maxTrim) {
          note = `試行時間が短いため、前後除外を ${trimSec.toFixed(2)}→${maxTrim.toFixed(2)} 秒に自動調整しました。`;
          trimSec = maxTrim;
        }
      }

      const t0 = time[0] + trimSec;
      const t1 = time[time.length - 1] - trimSec;
      let trimmed = trimByTime(time, copx, copy, t0, t1);

      if (trimmed.time.length < 20) {
        showMsg("前後除外後にデータが少なすぎます。前後除外(秒)を減らしてください。（欠損補完はセンタリング前に実施しています）", "warn");
        return;
      }

      // Centering (on raw COP)
      const cx = center(trimmed.copx);
      const cy = center(trimmed.copy);

      // Lowpass (filtfilt-style)
      const fx = filtfiltButter4Lowpass(cx, fs, cutoff);
      const fy = filtfiltButter4Lowpass(cy, fs, cutoff);

      // Unit conversion after filter
      const xcm = fx.map(v => v * unitToCm);
      const ycm = fy.map(v => v * unitToCm);

      // 95% ellipse area
      const cov = cov2(xcm, ycm);
      const ell = ellipseFromCov(cov);
      const ellipseArea = ell.areaCm2;

      // Convex hull area
      const hull = convexHull(xcm, ycm);
      const hullArea = polygonArea(hull);

      // Welch PSD & band ratios
      const bands = {
        LF: [0.0, 0.3],
        MF: [0.3, 1.0],
        HF: [1.0, 3.0],
      };

      const psdX = welchPSD(fx, fs);
      const psdY = welchPSD(fy, fs);

      // combined power density
      const fgrid = psdX.f;
      const pTotal = fgrid.map((_, i) => psdX.p[i] + psdY.p[i]);

      const ratiosTotal = bandRatiosFromPSD(fgrid, pTotal, bands);
      const ratiosML = bandRatiosFromPSD(psdX.f, psdX.p, bands);
      const ratiosAP = bandRatiosFromPSD(psdY.f, psdY.p, bands);

      // "MF" and "RPF" – users often want MF(0.3–1.0) and HF(1–3) ratio.
      // We show all three, but highlight MF and HF (as RPF).
      const MF = ratiosTotal.MF;
      const RPF = ratiosTotal.HF; // displayed as RPF(1–3Hz)

      const res = {
        fs, cutoff, trimSec, unit, unitToCm,
        autoTrim, resample,
        note,
        n: trimmed.time.length,
        duration: trimmed.time[trimmed.time.length - 1] - trimmed.time[0],
        ellipseAreaCm2: ellipseArea,
        hullAreaCm2: hullArea,
        ratiosTotal,
        ratiosML,
        ratiosAP,
        ellipse: ell, // points for display
        hull,
        xcm,
        ycm,
        fgrid,
      };
      lastResult = res;

      showMsg("解析完了。", "ok");
      renderMetrics(res, MF, RPF);
      await renderPlots(res, MF, RPF);

      dlMetricsBtn.disabled = false;
      dlBandsBtn.disabled = false;
      dlPngBtn.disabled = false;

      dlMetricsBtn.onclick = () => downloadText(buildMetricsCsv(res, MF, RPF), `${lastFilename}_metrics.csv`, "text/csv");
      dlBandsBtn.onclick = () => downloadText(buildBandsCsv(res), `${lastFilename}_bands.csv`, "text/csv");
      dlPngBtn.onclick = async () => {
        const node = $("plotCop");
        const imgData = await Plotly.toImage(node, {format:"png", height: 700, width: 900, scale: 2});
        downloadDataUrl(imgData, `${lastFilename}_cop.png`);
      };

    } catch (e) {
      console.error(e);
      showMsg(`エラー：${e?.message || e}`, "warn");
    }
  });

  function resetOutputs() {
    lastResult = null;
    msgEl.style.display = "none";
    metricsEl.style.display = "none";
    metricsEl.innerHTML = "";
    detailsEl.textContent = "";
    dlMetricsBtn.disabled = true;
    dlBandsBtn.disabled = true;
    dlPngBtn.disabled = true;
    Plotly.purge("plotCop");
    Plotly.purge("plotBands");
  }

  function showMsg(text, kind) {
    msgEl.textContent = text;
    msgEl.className = (kind === "ok") ? "ok" : "warn";
    msgEl.style.display = "block";
  }

  // --- CSV parsing ---
  function parseCsvFile(file) {
    return new Promise((resolve, reject) => {
      Papa.parse(file, {
        header: true,
        skipEmptyLines: true,
        dynamicTyping: false,
        complete: (results) => {
          if (results.errors && results.errors.length) {
            reject(new Error(results.errors[0].message || "CSV parse error"));
          } else {
            resolve(results.data || []);
          }
        },
      });
    });
  }

  
  function toNumber(v) {
    // Robust parsing: handles trailing units, commas as decimal separators, etc.
    if (v === null || v === undefined) return NaN;
    const s0 = String(v).trim();
    if (!s0) return NaN;

    // If looks like "1,23" (comma decimal) and no dot, convert comma to dot.
    let s = s0;
    if (s.includes(",") && !s.includes(".")) {
      // avoid treating thousands separators in cases like "1,234"
      const parts = s.split(",");
      if (parts.length === 2 && parts[1].length <= 3) s = parts[0] + "." + parts[1];
    }

    // Extract first float-like token
    const m = s.match(/[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?/);
    if (!m) return NaN;
    const num = Number(m[0]);
    return Number.isFinite(num) ? num : NaN;
  }

  function medianDiff(arr) {
    if (arr.length < 3) return NaN;
    const diffs = [];
    for (let i = 1; i < arr.length; i++) {
      const d = arr[i] - arr[i-1];
      if (Number.isFinite(d) && d > 0) diffs.push(d);
    }
    if (!diffs.length) return NaN;
    diffs.sort((a,b)=>a-b);
    return diffs[Math.floor(diffs.length/2)];
  }

  function normalizeTimeUnits(time, fs) {
    // Heuristics: if time is in ms or sample index, convert to seconds.
    // Returns { timeSec, note }
    const n = time.length;
    if (n < 2) return { timeSec: time.slice(), note: "" };

    const t0 = time[0], t1 = time[n-1];
    const durObs = t1 - t0;
    if (!Number.isFinite(durObs) || durObs <= 0) return { timeSec: time.slice(), note: "" };

    const durExp = (n - 1) / fs; // expected seconds if time is already seconds at fs
    if (!Number.isFinite(durExp) || durExp <= 0) return { timeSec: time.slice(), note: "" };

    const ratio = durObs / durExp;
    const dtMed = medianDiff(time);

    // milliseconds: ratio ~ 1000
    if (ratio > 700 && ratio < 1300) {
      return { timeSec: time.map(t => t / 1000), note: "time列をms→秒に自動変換しました。" };
    }

    // sample index: dt ~ 1 and ratio ~ fs
    if (Number.isFinite(dtMed) && dtMed > 0.8 && dtMed < 1.2 && ratio > 0.7*fs && ratio < 1.3*fs) {
      return { timeSec: time.map(t => t / fs), note: "time列をサンプル番号とみなし、fsで秒に換算しました。" };
    }

    // already seconds
    return { timeSec: time.slice(), note: "" };
  }

function normalizeKey(k) {
    return String(k || "").trim().toLowerCase().replace(/\s+/g, "");
  }

  function extractColumns(rows) {
    if (!rows.length) return { ok: false, error: "CSVが空です。" };

    // Map headers
    const headers = Object.keys(rows[0] || {});
    const normToOrig = new Map();
    for (const h of headers) normToOrig.set(normalizeKey(h), h);

    const pick = (...cands) => {
      for (const c of cands) {
        const k = normalizeKey(c);
        if (normToOrig.has(k)) return normToOrig.get(k);
      }
      return null;
    };

    const timeKey = pick("time", "t", "timestamp", "sec", "seconds");
    const xKey = pick("copx", "cop_x", "cop-x", "ml", "x");
    const yKey = pick("copy", "copyy", "cop_y", "cop-y", "ap", "y");

    if (!timeKey || !xKey || !yKey) {
      return { ok: false, error: `必須列が見つかりません。見つかった列：${headers.join(", ")}` };
    }

    const time = [];
    const copx = [];
    const copy = [];

    for (const r of rows) {
      const t = toNumber(r[timeKey]);
      const x = toNumber(r[xKey]);
      const y = toNumber(r[yKey]);
      if (Number.isFinite(t) && Number.isFinite(x) && Number.isFinite(y)) {
        time.push(t); copx.push(x); copy.push(y);
      }
    }
    if (time.length < 2) return { ok: false, error: "数値として読める行がありません。" };

    return { ok: true, time, copx, copy, keys: { timeKey, xKey, yKey } };
  }

  
  // --- preprocessing

  function sortByTime(time, x, y) {
    const idx = time.map((t,i)=>i).sort((a,b)=>time[a]-time[b]);
    const nt=[], nx=[], ny=[];
    let lastT = null;
    for (const i of idx) {
      const t = time[i];
      if (!Number.isFinite(t)) continue;
      // de-duplicate exact same timestamp: keep the last
      if (lastT !== null && t === lastT) {
        nt[nt.length-1] = t;
        nx[nx.length-1] = x[i];
        ny[ny.length-1] = y[i];
      } else {
        nt.push(t); nx.push(x[i]); ny.push(y[i]);
        lastT = t;
      }
    }
    return { time: nt, copx: nx, copy: ny };
  }

  function resampleLinearUniform(time, x, y, fs) {
    if (time.length < 2) return { time: time.slice(), copx: x.slice(), copy: y.slice() };
    const tStart = time[0];
    const tEnd = time[time.length - 1];
    const dt = 1 / fs;
    const n = Math.max(2, Math.floor((tEnd - tStart) / dt) + 1);
    const rt = new Array(n);
    const rx = new Array(n);
    const ry = new Array(n);
    let j = 0;

    for (let i = 0; i < n; i++) {
      const t = tStart + i * dt;
      rt[i] = t;

      while (j < time.length - 2 && time[j + 1] < t) j++;

      const tA = time[j];
      const tB = time[j + 1];
      const xA = x[j], xB = x[j + 1];
      const yA = y[j], yB = y[j + 1];

      if (tB === tA) {
        rx[i] = xA;
        ry[i] = yA;
      } else {
        const w = (t - tA) / (tB - tA);
        const ww = Math.max(0, Math.min(1, w));
        rx[i] = xA + ww * (xB - xA);
        ry[i] = yA + ww * (yB - yA);
      }
    }
    return { time: rt, copx: rx, copy: ry };
  }

  function trimByTime(time, x, y, t0, t1) {
    const nt = [];
    const nx = [];
    const ny = [];
    for (let i = 0; i < time.length; i++) {
      const t = time[i];
      if (t >= t0 && t <= t1) {
        nt.push(t); nx.push(x[i]); ny.push(y[i]);
      }
    }
    return { time: nt, copx: nx, copy: ny };
  }

  function mean(arr) {
    let s = 0;
    for (const v of arr) s += v;
    return s / arr.length;
  }

  function center(arr) {
    const m = mean(arr);
    return arr.map(v => v - m);
  }

  // --- Butterworth 4th-order lowpass at Wn=0.2 (cutoff=10Hz, fs=100Hz)
  // To closely match the notebook, we hardcode SciPy's butter(4, 0.2) coefficients,
  // and implement filtfilt-style forward/backward filtering with odd extension.
  //
  // If user changes fs/cutoff, we still apply these coefficients as long as the ratio is 0.2.
  // For other ratios, we fall back to a simple biquad cascade approximation (still zero-phase).
  const B_FIXED = [0.004824343357716228, 0.01929737343086491, 0.028946060146297367, 0.01929737343086491, 0.004824343357716228];
  const A_FIXED = [1.0, -2.369513007182037, 2.3139884144158806, -1.0546654058785673, 0.18737949236818535];
  const ZI_FIXED = [0.9951756566422838, -1.3936347201975356, 0.891407628579869, -0.1825551467889937];

  function approxEqual(a, b, eps=1e-6) { return Math.abs(a-b) <= eps; }

  function filtfiltButter4Lowpass(x, fs, cutoff) {
    const n = x.length;
    if (n <= 15) return x.slice();

    const wn = cutoff / (0.5 * fs);
    // If exactly the notebook condition (wn=0.2), use fixed coefficients
    if (approxEqual(wn, 0.2, 1e-3)) {
      return filtfiltIIR(x, B_FIXED, A_FIXED, ZI_FIXED);
    }

    // Otherwise: use a bilinear-transform based RBJ biquad cascade (2x biquad) to approximate 4th-order.
    // This is a pragmatic fallback; default UI values match wn=0.2.
    const biquads = designButterworth4LowpassBiquads(fs, cutoff);
    let y = x.slice();
    for (const bq of biquads) y = filtfiltBiquad(y, bq);
    return y;
  }

  function filtfiltIIR(x, b, a, ziBase) {
    const padlen = 3 * (Math.max(a.length, b.length) - 1); // 12
    if (x.length <= padlen) return x.slice();

    const ext = oddExtend(x, padlen);

    // forward
    const y1 = lfilter(ext, b, a, ziBase, ext[0]);

    // backward
    const rev = y1.slice().reverse();
    const y2 = lfilter(rev, b, a, ziBase, rev[0]).reverse();

    // trim
    return y2.slice(padlen, padlen + x.length);
  }

  function oddExtend(x, padlen) {
    // ext = [2*x0 - x[1:padlen+1][::-1], x, 2*xN - x[-padlen-1:-1][::-1]]
    const n = x.length;
    const left = [];
    for (let i = 1; i <= padlen; i++) left.push(2*x[0] - x[i]);
    left.reverse();

    const right = [];
    for (let i = 1; i <= padlen; i++) right.push(2*x[n-1] - x[n-1-i]);
    right.reverse();

    return left.concat(x, right);
  }

  function lfilter(x, b, a, ziBase, x0) {
    // Direct Form II Transposed with initial state scaled by x0 (steady-state for step)
    const nb = b.length;
    const na = a.length;
    const order = Math.max(nb, na) - 1;
    const z = new Array(order).fill(0);
    for (let i = 0; i < order; i++) z[i] = ziBase[i] * x0;

    const y = new Array(x.length).fill(0);
    for (let n = 0; n < x.length; n++) {
      const xn = x[n];
      let yn = b[0]*xn + z[0];
      for (let i = 1; i < order; i++) {
        z[i-1] = (b[i] || 0)*xn + z[i] - (a[i] || 0)*yn;
      }
      z[order-1] = (b[order] || 0)*xn - (a[order] || 0)*yn;
      y[n] = yn;
    }
    return y;
  }

  // --- Biquad fallback (RBJ) ---
  function designButterworth4LowpassBiquads(fs, fc) {
    // Two cascaded 2nd-order Butterworth sections with Q values:
    // For 4th order: Q1=0.5411961, Q2=1.30656296 (from analog prototype)
    const Qs = [0.5411961, 1.3065630];
    return Qs.map(Q => designRBJLowpass(fs, fc, Q));
  }

  function designRBJLowpass(fs, fc, Q) {
    const w0 = 2*Math.PI*fc/fs;
    const cosw0 = Math.cos(w0);
    const sinw0 = Math.sin(w0);
    const alpha = sinw0/(2*Q);

    let b0 = (1 - cosw0)/2;
    let b1 = 1 - cosw0;
    let b2 = (1 - cosw0)/2;
    let a0 = 1 + alpha;
    let a1 = -2*cosw0;
    let a2 = 1 - alpha;

    // normalize
    b0/=a0; b1/=a0; b2/=a0; a1/=a0; a2/=a0;
    return { b0,b1,b2,a1,a2 };
  }

  function filtfiltBiquad(x, bq) {
    const padlen = 12;
    if (x.length <= padlen) return x.slice();
    const ext = oddExtend(x, padlen);
    const y1 = biquadFilter(ext, bq);
    const y2 = biquadFilter(y1.slice().reverse(), bq).reverse();
    return y2.slice(padlen, padlen + x.length);
  }

  function biquadFilter(x, bq) {
    const y = new Array(x.length);
    let z1=0, z2=0;
    for (let i=0;i<x.length;i++) {
      const xn = x[i];
      const yn = bq.b0*xn + z1;
      z1 = bq.b1*xn + z2 - bq.a1*yn;
      z2 = bq.b2*xn - bq.a2*yn;
      y[i]=yn;
    }
    return y;
  }

  // --- covariance, ellipse ---
  function cov2(x, y) {
    const n = x.length;
    const mx = mean(x), my = mean(y);
    let sxx=0, syy=0, sxy=0;
    for (let i=0;i<n;i++) {
      const dx = x[i]-mx;
      const dy = y[i]-my;
      sxx += dx*dx;
      syy += dy*dy;
      sxy += dx*dy;
    }
    const denom = (n-1);
    return [[sxx/denom, sxy/denom],[sxy/denom, syy/denom]];
  }

  function ellipseFromCov(cov) {
    // cov = [[a,b],[b,d]]
    const a = cov[0][0], b = cov[0][1], d = cov[1][1];
    const tr = a + d;
    const disc = Math.sqrt((a - d)*(a - d) + 4*b*b);
    const l1 = (tr + disc)/2;
    const l2 = (tr - disc)/2;

    // eigenvector for l1
    let vx, vy;
    if (Math.abs(b) > 1e-12) {
      vx = l1 - d;
      vy = b;
    } else {
      // diagonal
      vx = 1; vy = 0;
      if (a < d) { vx = 0; vy = 1; }
    }
    const norm = Math.hypot(vx, vy);
    vx/=norm; vy/=norm;
    const angle = Math.atan2(vy, vx);

    // chi2.ppf(0.95,2) = 5.991464547107979
    const chi2ppf = 5.991464547107979;
    const scale = Math.sqrt(chi2ppf);

    const A = scale * Math.sqrt(Math.max(l1, 0));
    const B = scale * Math.sqrt(Math.max(l2, 0));

    const areaCm2 = Math.PI * chi2ppf * Math.sqrt(Math.max(l1,0) * Math.max(l2,0));

    // ellipse points (center at 0,0)
    const pts = [];
    const N=400;
    for (let i=0;i<=N;i++) {
      const t = (2*Math.PI*i)/N;
      const ex = A*Math.cos(t);
      const ey = B*Math.sin(t);
      // rotate
      const rx = ex*Math.cos(angle) - ey*Math.sin(angle);
      const ry = ex*Math.sin(angle) + ey*Math.cos(angle);
      pts.push([rx, ry]);
    }
    return { areaCm2, angleRad: angle, axisA: A, axisB: B, points: pts };
  }

  // --- convex hull (Monotonic chain) ---
  function convexHull(x, y) {
    const pts = x.map((vx, i) => [vx, y[i]]);
    // Sort by x then y
    pts.sort((p,q) => (p[0]-q[0]) || (p[1]-q[1]));
    const cross = (o,a,b) => (a[0]-o[0])*(b[1]-o[1]) - (a[1]-o[1])*(b[0]-o[0]);

    const lower = [];
    for (const p of pts) {
      while (lower.length >= 2 && cross(lower[lower.length-2], lower[lower.length-1], p) <= 0) {
        lower.pop();
      }
      lower.push(p);
    }
    const upper = [];
    for (let i = pts.length - 1; i >= 0; i--) {
      const p = pts[i];
      while (upper.length >= 2 && cross(upper[upper.length-2], upper[upper.length-1], p) <= 0) {
        upper.pop();
      }
      upper.push(p);
    }
    upper.pop();
    lower.pop();
    return lower.concat(upper); // closed polygon not repeated
  }

  function polygonArea(poly) {
    if (!poly || poly.length < 3) return 0;
    let s=0;
    for (let i=0;i<poly.length;i++) {
      const [x1,y1]=poly[i];
      const [x2,y2]=poly[(i+1)%poly.length];
      s += x1*y2 - x2*y1;
    }
    return Math.abs(s)/2;
  }

  // --- Welch PSD (Hann, 50% overlap) ---
  function hann(N) {
    const w = new Array(N);
    for (let n=0;n<N;n++) w[n] = 0.5 - 0.5*Math.cos(2*Math.PI*n/(N-1));
    return w;
  }

  function welchPSD(x, fs) {
    const N = x.length;
    const nperseg = Math.min(1024, N);
    const noverlap = Math.floor(nperseg/2);
    const step = nperseg - noverlap;
    if (step <= 0) throw new Error("noverlap too large");

    const w = hann(nperseg);
    let w2sum = 0;
    for (const v of w) w2sum += v*v;

    const fft = new FFT(nperseg);
    const out = fft.createComplexArray();
    const input = fft.createComplexArray();

    const nfreq = Math.floor(nperseg/2) + 1;
    const P = new Array(nfreq).fill(0);
    let nseg = 0;

    for (let start=0; start + nperseg <= N; start += step) {
      // segment
      let m=0;
      for (let i=0;i<nperseg;i++) m += x[start+i];
      m /= nperseg;

      for (let i=0;i<nperseg;i++) {
        const v = (x[start+i] - m) * w[i];
        input[2*i] = v;
        input[2*i + 1] = 0;
      }

      fft.transform(out, input);

      // one-sided power density
      // scale = 1 / (fs * sum(w^2))
      const scale = 1 / (fs * w2sum);

      for (let k=0;k<nfreq;k++) {
        const re = out[2*k];
        const im = out[2*k + 1];
        let pk = (re*re + im*im) * scale;

        // double non-DC and non-Nyquist
        if (k !== 0 && k !== nperseg/2) pk *= 2;

        P[k] += pk;
      }
      nseg += 1;
    }

    if (nseg === 0) throw new Error("データが短すぎてWelch計算できません。");

    for (let k=0;k<P.length;k++) P[k] /= nseg;

    const f = new Array(nfreq);
    for (let k=0;k<nfreq;k++) f[k] = (k * fs) / nperseg;
    return { f, p: P, nperseg, noverlap, nseg };
  }

  function trapz(x, y) {
    // trapezoidal integral over y(x)
    let s = 0;
    for (let i=0;i<x.length-1;i++) {
      const dx = x[i+1]-x[i];
      s += dx * (y[i] + y[i+1]) / 2;
    }
    return s;
  }

  function bandRatiosFromPSD(f, p, bands) {
    // total power within 0–3 Hz
    const fmin = 0.0, fmax = 3.0;
    const idx = [];
    for (let i=0;i<f.length;i++) if (f[i] >= fmin && f[i] <= fmax) idx.push(i);
    if (idx.length < 2) return { LF: NaN, MF: NaN, HF: NaN };

    const ft = idx.map(i => f[i]);
    const pt = idx.map(i => p[i]);
    const total = trapz(ft, pt);
    const ratios = {};
    for (const [name, [a,b]] of Object.entries(bands)) {
      const j = [];
      for (let i=0;i<f.length;i++) if (f[i] >= a && f[i] <= b) j.push(i);
      if (j.length < 2 || !(total > 0)) {
        ratios[name] = NaN;
      } else {
        const fb = j.map(i => f[i]);
        const pb = j.map(i => p[i]);
        ratios[name] = trapz(fb, pb) / total * 100;
      }
    }
    return ratios;
  }

  // --- rendering ---
  function fmt(v, digits=2) {
    if (!Number.isFinite(v)) return "—";
    return v.toFixed(digits);
  }

  function renderMetrics(res, MF, RPF) {
    const items = [
      { k: "楕円面積 (95%)", v: fmt(res.ellipseAreaCm2, 3), u: "cm²" },
      { k: "凸包面積", v: fmt(res.hullAreaCm2, 3), u: "cm²" },
      { k: "MF (0.3–1.0Hz)", v: fmt(MF, 2), u: "%" },
      { k: "RPF (1–3Hz)", v: fmt(RPF, 2), u: "%" },
      { k: "LF (0–0.3Hz)", v: fmt(res.ratiosTotal.LF, 2), u: "%" },
      { k: "データ点数", v: String(res.n), u: "" },
    ];

    metricsEl.innerHTML = items.map(it => `
      <div class="metric">
        <div class="k">${it.k}</div>
        <div class="v">${it.v}<span class="u">${it.u}</span></div>
      </div>
    `).join("");
    metricsEl.style.display = "grid";

    detailsEl.textContent =
      `設定: FS=${res.fs}Hz / LPF=${res.cutoff}Hz / 前後除外=${res.trimSec}s / 単位=${res.unit} / 補完=${res.resample ? "ON" : "OFF"} / Welch(nperseg=${Math.min(1024,res.n)})` + (res.note ? `
${res.note}` : "");
  }

  async function renderPlots(res, MF, RPF) {
    // --- COP plot with ellipse + hull ---
    const hullX = res.hull.map(p => p[0]).concat(res.hull[0]?.[0] ?? []);
    const hullY = res.hull.map(p => p[1]).concat(res.hull[0]?.[1] ?? []);

    const ellX = res.ellipse.points.map(p => p[0]);
    const ellY = res.ellipse.points.map(p => p[1]);

    const traceCop = {
      x: res.xcm, y: res.ycm, mode: "lines", name: "COP",
      line: { width: 2 },
      opacity: 0.8,
    };
    const traceHull = {
      x: hullX, y: hullY, mode: "lines", name: "凸包",
      line: { width: 2, dash: "solid" },
      opacity: 0.9,
    };
    const traceEllipse = {
      x: ellX, y: ellY, mode: "lines", name: "95%楕円",
      line: { width: 3, dash: "dash" },
      opacity: 0.9,
    };
    const traceCenter = {
      x: [0], y: [0], mode:"markers", name:"中心(0,0)",
      marker: { size: 10, symbol: "circle-open", line: { width: 2 } }
    };

    const layoutCop = {
      title: { text: `COP（楕円面積=${fmt(res.ellipseAreaCm2,2)} cm² / 凸包面積=${fmt(res.hullAreaCm2,2)} cm²）`, font: { size: 14 } },
      xaxis: { title: "ML (cm)", zeroline: true, scaleanchor: "y", scaleratio: 1 },
      yaxis: { title: "AP (cm)", zeroline: true },
      margin: { l: 50, r: 20, t: 60, b: 50 },
      legend: { orientation: "h" },
    };

    await Plotly.newPlot("plotCop", [traceCop, traceHull, traceEllipse, traceCenter], layoutCop, {responsive:true});

    // --- Band ratios bar plot (ML/AP + Total) ---
    const labels = ["LF", "MF", "HF"];
    const totalVals = labels.map(k => res.ratiosTotal[k]);
    const mlVals = labels.map(k => res.ratiosML[k]);
    const apVals = labels.map(k => res.ratiosAP[k]);

    const tTotal = { x: labels, y: totalVals, type:"bar", name:"Total(ML+AP)" };
    const tML = { x: labels, y: mlVals, type:"bar", name:"ML" };
    const tAP = { x: labels, y: apVals, type:"bar", name:"AP" };

    const layoutBands = {
      title: { text: `帯域比（0–3Hz=100%）：MF=${fmt(MF,2)}% / RPF(=HF)=${fmt(RPF,2)}%`, font: { size: 14 } },
      yaxis: { title: "Power Ratio (%)", range: [0, 100] },
      barmode: "group",
      margin: { l: 55, r: 20, t: 60, b: 40 },
      legend: { orientation:"h" },
    };

    await Plotly.newPlot("plotBands", [tTotal, tML, tAP], layoutBands, {responsive:true});
  }

  // --- downloads ---
  function downloadText(text, filename, mime) {
    const blob = new Blob([text], { type: mime || "text/plain" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = filename;
    a.click();
    URL.revokeObjectURL(url);
  }

  function downloadDataUrl(dataUrl, filename) {
    const a = document.createElement("a");
    a.href = dataUrl;
    a.download = filename;
    a.click();
  }

  function buildMetricsCsv(res, MF, RPF) {
    const lines = [];
    lines.push("key,value,unit");
    lines.push(`ellipse_area_95,${res.ellipseAreaCm2},cm2`);
    lines.push(`convex_hull_area,${res.hullAreaCm2},cm2`);
    lines.push(`LF_ratio,${res.ratiosTotal.LF},percent`);
    lines.push(`MF_ratio,${MF},percent`);
    lines.push(`RPF_ratio_HF_1_3,${RPF},percent`);
    lines.push(`HF_ratio,${res.ratiosTotal.HF},percent`);
    lines.push(`n_points,${res.n},count`);
    lines.push(`fs,${res.fs},Hz`);
    lines.push(`lowpass,${res.cutoff},Hz`);
    lines.push(`trim_sec,${res.trimSec},sec`);
    lines.push(`unit,${res.unit},`);
    return lines.join("\n");
  }

  function buildBandsCsv(res) {
    const lines = [];
    lines.push("band,total_percent,ML_percent,AP_percent");
    for (const band of ["LF","MF","HF"]) {
      lines.push(`${band},${res.ratiosTotal[band]},${res.ratiosML[band]},${res.ratiosAP[band]}`);
    }
    return lines.join("\n");
  }
})();
