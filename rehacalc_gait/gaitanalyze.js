/* ========================================================
   Gait Acceleration Analysis  –  gaitanalyze.js
   L3 腰椎加速度データから歩行指標を算出
   Based on R code: ACC_data_analysis.r
   ======================================================== */

/* ---------- globals ---------- */
let rawData = null;

/* ========== UI wiring ========== */
const $  = id => document.getElementById(id);
const fs = () => parseInt($('fs').value, 10) || 100;

$('file').addEventListener('change', () => { $('run').disabled = !$('file').files.length; });
$('run').addEventListener('click', run);

/* ========== main pipeline ========== */
async function run() {
  $('run').disabled = true;
  $('msg').style.display = 'none';
  $('metrics').style.display = 'none';
  $('details').textContent = '';

  try {
    const text = await $('file').files[0].text();
    const parsed = Papa.parse(text.trim(), { header: true, dynamicTyping: true, skipEmptyLines: true });
    if (parsed.errors.length) throw new Error('CSV parse error: ' + parsed.errors[0].message);

    const rows = parsed.data;
    if (rows.length < 10) throw new Error('データが少なすぎます');

    /* --- detect columns --- */
    const cols = detectColumns(parsed.meta.fields, rows);

    /* --- extract arrays --- */
    let time = rows.map(r => r[cols.time]);
    let ap   = rows.map(r => r[cols.ap]);
    let ml   = rows.map(r => r[cols.ml]);
    let vt   = rows.map(r => r[cols.vt]);

    /* remove NaN */
    const valid = time.map((t, i) => !isNaN(t) && !isNaN(ap[i]) && !isNaN(ml[i]) && !isNaN(vt[i]));
    time = time.filter((_, i) => valid[i]);
    ap   = ap.filter((_, i) => valid[i]);
    ml   = ml.filter((_, i) => valid[i]);
    vt   = vt.filter((_, i) => valid[i]);

    /* --- auto-detect sampling rate from time column --- */
    const userFs = parseInt($('fs').value, 10);
    let sampleRate;
    if ($('autoFs').checked) {
      const diffs = [];
      const nSample = Math.min(time.length - 1, 500);
      for (let i = 1; i <= nSample; i++) diffs.push(time[i] - time[i - 1]);
      diffs.sort((a, b) => a - b);
      const medianDiff = diffs[Math.floor(diffs.length / 2)];
      sampleRate = Math.round(1 / medianDiff);
      $('fs').value = sampleRate;
    } else {
      sampleRate = userFs || 100;
    }

    /* --- resample to uniform intervals if needed --- */
    const dt = 1 / sampleRate;
    const tStart = time[0];
    const tEnd = time[time.length - 1];
    const nUniform = Math.floor((tEnd - tStart) / dt) + 1;
    if (nUniform > 10) {
      const uTime = Array.from({ length: nUniform }, (_, i) => tStart + i * dt);
      ap = linearInterp(time, ap, uTime);
      ml = linearInterp(time, ml, uTime);
      vt = linearInterp(time, vt, uTime);
      time = uTime;
    }

    /* --- 1) Butterworth low-pass filter 17 Hz --- */
    const cutoff = parseFloat($('cutoff').value) || 17;
    const coeffs = butterworth4(cutoff, sampleRate);
    const filt_ap = filtfilt(coeffs, ap);
    const filt_ml = filtfilt(coeffs, ml);
    const filt_vt = filtfilt(coeffs, vt);

    /* --- 2) Walking bout detection (moving SD) --- */
    const winSize = Math.round(0.3 * sampleRate);
    const msd_ap = movingSD(filt_ap, winSize);
    const msd_ml = movingSD(filt_ml, winSize);
    const msd_vt = movingSD(filt_vt, winSize);
    const msd_sum = msd_ap.map((v, i) => v + msd_ml[i] + msd_vt[i]);

    const threshold = parseFloat($('threshold').value) || 0.03;
    const bout = detectBout(msd_sum, threshold);

    if (bout.end - bout.start < sampleRate) throw new Error('歩行区間が短すぎます');

    const bout_ap = filt_ap.slice(bout.start, bout.end + 1);
    const bout_ml = filt_ml.slice(bout.start, bout.end + 1);
    const bout_vt = filt_vt.slice(bout.start, bout.end + 1);
    const boutDuration = bout_ap.length / sampleRate;

    /* --- 3) Stride detection --- */
    const coeffsLow = butterworth4(2, sampleRate);
    const filt_bout_ap = filtfilt(coeffsLow, bout_ap);
    const peaks = findPeaks(filt_bout_ap, Math.round(sampleRate / 4));

    const exceptStride = parseInt($('exceptStride').value, 10) || 2;

    if (peaks.length < 2 * exceptStride + 3) throw new Error('ストライドが検出できません（ピーク不足）');

    /* --- 4) Per-stride metrics --- */
    const strideTimes = [];
    const HR_ap_arr = [], HR_ml_arr = [], HR_vt_arr = [];
    const rms_ap_arr = [], rms_ml_arr = [];

    const nStrides = Math.floor(peaks.length / 2);
    for (let p = exceptStride; p < nStrides - exceptStride; p++) {
      const i0 = peaks[2 * p - 1]; // using 0-based peak indices
      const i1 = peaks[2 * p + 1];
      if (i0 === undefined || i1 === undefined) continue;

      const cyc_ap = bout_ap.slice(i0, i1 + 1);
      const cyc_ml = bout_ml.slice(i0, i1 + 1);
      const cyc_vt = bout_vt.slice(i0, i1 + 1);

      if (cyc_ap.length < 4) continue;

      strideTimes.push(cyc_ap.length / sampleRate);

      /* Harmonic Ratio via FFT */
      const hr = calcHarmonicRatios(cyc_ap, cyc_ml, cyc_vt, sampleRate);
      HR_ap_arr.push(hr.ap);
      HR_ml_arr.push(hr.ml);
      HR_vt_arr.push(hr.vt);

      /* RMS of each cycle */
      rms_ap_arr.push(rms(cyc_ap));
      rms_ml_arr.push(rms(cyc_ml));
    }

    if (strideTimes.length === 0) throw new Error('有効なストライドが検出されませんでした');

    /* --- 5) Walking speed estimate --- */
    /* 10m歩行距離を入力させるか、CSVのtime列から推定 */
    const walkDist = parseFloat($('walkDist').value) || 10; // metres
    const walkSpeed = walkDist / boutDuration; // m/s

    /* --- 6) Aggregate metrics --- */
    const strideMean   = mean(strideTimes);
    const strideSD     = sd(strideTimes);
    const strideCV     = strideSD / strideMean * 100;
    const cadence      = 60 / strideMean * 2; // steps/min (2 steps per stride)
    const hrAP         = mean(HR_ap_arr);
    const hrML         = mean(HR_ml_arr);
    const hrVT         = mean(HR_vt_arr);
    const rmsAP        = mean(rms_ap_arr);
    const rmsML        = mean(rms_ml_arr);
    /* 歩行速度で正規化 */
    const rmsAP_norm   = rmsAP / walkSpeed;
    const rmsML_norm   = rmsML / walkSpeed;
    const trunkInstab  = Math.sqrt(rmsAP_norm ** 2 + rmsML_norm ** 2);

    /* --- 7) Display --- */
    showMetrics({
      walkSpeed, boutDuration, strideMean, strideSD, strideCV, cadence,
      hrAP, hrML, hrVT,
      rmsAP, rmsML, rmsAP_norm, rmsML_norm, trunkInstab,
      nStrides: strideTimes.length
    });

    /* --- 8) Plots --- */
    plotResults({
      time, filt_ap, filt_ml, filt_vt,
      bout, msd_sum, threshold,
      bout_ap, bout_ml, bout_vt,
      filt_bout_ap, peaks, exceptStride,
      sampleRate,
      strideTimes, HR_ap_arr, HR_ml_arr, HR_vt_arr
    });

    showMsg('解析完了', 'ok');
    $('dlCsv').disabled = false;
    $('dlPng').disabled = false;

    /* store for CSV export */
    rawData = {
      walkSpeed, boutDuration, strideMean, strideSD, strideCV, cadence,
      hrAP, hrML, hrVT, rmsAP, rmsML, rmsAP_norm, rmsML_norm, trunkInstab,
      nStrides: strideTimes.length,
      strideTimes, HR_ap_arr, HR_ml_arr, HR_vt_arr, rms_ap_arr, rms_ml_arr
    };

  } catch (e) {
    showMsg(e.message, 'warn');
    console.error(e);
  }
  $('run').disabled = false;
}

/* ========== Column auto-detection ========== */
function detectColumns(fields, rows) {
  const lower = fields.map(f => f.toLowerCase().trim());

  /* time column */
  let timeIdx = lower.findIndex(f => /^time/.test(f));
  if (timeIdx < 0) timeIdx = 0;

  /* Try axis assignment from user dropdowns first */
  const apSel = $('colAP').value;
  const mlSel = $('colML').value;
  const vtSel = $('colVT').value;
  if (apSel && mlSel && vtSel) {
    return { time: fields[timeIdx], ap: apSel, ml: mlSel, vt: vtSel };
  }

  /* Auto-detect: common patterns */
  let apIdx = -1, mlIdx = -1, vtIdx = -1;

  // L3 sensor convention: gFx=VT(上下), gFy=ML(左右), gFz=AP(前後)
  const hasGF = lower.some(f => /^gf[xyz]/.test(f));
  if (hasGF) {
    apIdx = lower.findIndex(f => /^gfz/.test(f));  // gFz = AP
    mlIdx = lower.findIndex(f => /^gfy/.test(f));  // gFy = ML
    vtIdx = lower.findIndex(f => /^gfx/.test(f));  // gFx = VT
  } else {
    apIdx = lower.findIndex(f => /ap|anteroposterior/.test(f));
    mlIdx = lower.findIndex(f => /ml|mediolateral/.test(f));
    vtIdx = lower.findIndex(f => /vt|vertical/.test(f));
  }

  // Fallback: columns 1,2,3 (after time)
  if (apIdx < 0) apIdx = timeIdx + 1;
  if (mlIdx < 0) mlIdx = timeIdx + 2;
  if (vtIdx < 0) vtIdx = timeIdx + 3;

  return {
    time: fields[timeIdx],
    ap: fields[Math.min(apIdx, fields.length - 1)],
    ml: fields[Math.min(mlIdx, fields.length - 1)],
    vt: fields[Math.min(vtIdx, fields.length - 1)]
  };
}

/* ========== Butterworth 4th order low-pass ========== */
function butterworth4(fc, fsSamp) {
  // Compute 4th-order Butterworth coefficients via bilinear transform
  const wc = Math.tan(Math.PI * fc / fsSamp);
  const wc2 = wc * wc;
  const wc3 = wc2 * wc;
  const wc4 = wc2 * wc2;

  // 4th order: cascade of two 2nd-order sections
  // Analog poles for 4th order Butterworth (damping coefficients, must be positive)
  const k1 = 2 * Math.cos(Math.PI * 3 / 8); // ≈ 0.7654
  const k2 = 2 * Math.cos(Math.PI * 1 / 8); // ≈ 1.8478

  // Section 1
  const a0_1 = 1 + k1 * wc + wc2;
  const b1 = [wc2 / a0_1, 2 * wc2 / a0_1, wc2 / a0_1];
  const a1 = [1, 2 * (wc2 - 1) / a0_1, (1 - k1 * wc + wc2) / a0_1];

  // Section 2
  const a0_2 = 1 + k2 * wc + wc2;
  const b2 = [wc2 / a0_2, 2 * wc2 / a0_2, wc2 / a0_2];
  const a2 = [1, 2 * (wc2 - 1) / a0_2, (1 - k2 * wc + wc2) / a0_2];

  // Combine into single 4th order (convolve numerator and denominator)
  const b = conv(b1, b2);
  const a = conv(a1, a2);

  return { b, a };
}

function conv(x, y) {
  const n = x.length + y.length - 1;
  const r = new Array(n).fill(0);
  for (let i = 0; i < x.length; i++)
    for (let j = 0; j < y.length; j++)
      r[i + j] += x[i] * y[j];
  return r;
}

/* ========== filtfilt (zero-phase filtering) ========== */
function filtfilt(coeffs, sig) {
  const { b, a } = coeffs;
  const padLen = 3 * Math.max(b.length, a.length);
  const n = sig.length;

  // Reflect-pad
  const padded = new Array(n + 2 * padLen);
  for (let i = 0; i < padLen; i++) {
    padded[i] = 2 * sig[0] - sig[padLen - i];
    padded[n + padLen + i] = 2 * sig[n - 1] - sig[n - 2 - i];
  }
  for (let i = 0; i < n; i++) padded[padLen + i] = sig[i];

  // Forward filter
  const fwd = iirFilter(b, a, padded);
  // Reverse
  fwd.reverse();
  const rev = iirFilter(b, a, fwd);
  rev.reverse();

  // Extract
  return rev.slice(padLen, padLen + n);
}

function iirFilter(b, a, x) {
  const n = x.length;
  const y = new Array(n).fill(0);
  for (let i = 0; i < n; i++) {
    y[i] = 0;
    for (let j = 0; j < b.length; j++) {
      if (i - j >= 0) y[i] += b[j] * x[i - j];
    }
    for (let j = 1; j < a.length; j++) {
      if (i - j >= 0) y[i] -= a[j] * y[i - j];
    }
    y[i] /= a[0];
  }
  return y;
}

/* ========== Moving standard deviation ========== */
function movingSD(arr, win) {
  const n = arr.length;
  const result = new Array(n).fill(0);
  const halfW = Math.floor(win / 2);
  for (let i = 0; i < n; i++) {
    const lo = Math.max(0, i - halfW);
    const hi = Math.min(n - 1, i + halfW);
    const seg = arr.slice(lo, hi + 1);
    const m = mean(seg);
    result[i] = Math.sqrt(seg.reduce((s, v) => s + (v - m) ** 2, 0) / seg.length);
  }
  return result;
}

/* ========== Walking bout detection ========== */
function detectBout(msdSum, thresh) {
  const binary = msdSum.map(v => v > thresh ? 1 : 0);

  // RLE to find longest run of 1s
  const runs = [];
  let cur = binary[0], count = 0;
  for (let i = 0; i < binary.length; i++) {
    if (binary[i] === cur) {
      count++;
    } else {
      runs.push({ val: cur, len: count });
      cur = binary[i];
      count = 1;
    }
  }
  runs.push({ val: cur, len: count });

  // Find longest run of 1
  let maxLen = 0, maxIdx = 0;
  for (let i = 0; i < runs.length; i++) {
    if (runs[i].val === 1 && runs[i].len > maxLen) {
      maxLen = runs[i].len;
      maxIdx = i;
    }
  }

  let start = 0;
  for (let i = 0; i < maxIdx; i++) start += runs[i].len;
  return { start, end: start + maxLen - 1 };
}

/* ========== Peak finding ========== */
function findPeaks(arr, minDist) {
  const peaks = [];
  for (let i = 1; i < arr.length - 1; i++) {
    if (arr[i] > arr[i - 1] && arr[i] >= arr[i + 1]) {
      peaks.push(i);
    }
  }
  // Enforce minimum distance: keep tallest peaks
  const filtered = [];
  for (const p of peaks) {
    if (filtered.length === 0 || p - filtered[filtered.length - 1] >= minDist) {
      filtered.push(p);
    } else {
      // Keep the taller one
      const lastIdx = filtered[filtered.length - 1];
      if (arr[p] > arr[lastIdx]) {
        filtered[filtered.length - 1] = p;
      }
    }
  }
  return filtered;
}

/* ========== FFT (Cooley-Tukey radix-2, with zero-padding) ========== */
function fft(re) {
  const N0 = re.length;
  // Zero-pad to next power of 2
  let N = 1;
  while (N < N0) N <<= 1;
  const real = new Float64Array(N);
  const imag = new Float64Array(N);
  for (let i = 0; i < N0; i++) real[i] = re[i];

  // Bit-reversal
  for (let i = 1, j = 0; i < N; i++) {
    let bit = N >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j) {
      [real[i], real[j]] = [real[j], real[i]];
      [imag[i], imag[j]] = [imag[j], imag[i]];
    }
  }

  for (let len = 2; len <= N; len <<= 1) {
    const ang = -2 * Math.PI / len;
    const wRe = Math.cos(ang), wIm = Math.sin(ang);
    for (let i = 0; i < N; i += len) {
      let curRe = 1, curIm = 0;
      for (let j = 0; j < len / 2; j++) {
        const uRe = real[i + j], uIm = imag[i + j];
        const vRe = real[i + j + len / 2] * curRe - imag[i + j + len / 2] * curIm;
        const vIm = real[i + j + len / 2] * curIm + imag[i + j + len / 2] * curRe;
        real[i + j] = uRe + vRe;
        imag[i + j] = uIm + vIm;
        real[i + j + len / 2] = uRe - vRe;
        imag[i + j + len / 2] = uIm - vIm;
        const newCurRe = curRe * wRe - curIm * wIm;
        curIm = curRe * wIm + curIm * wRe;
        curRe = newCurRe;
      }
    }
  }

  // Amplitude spectrum (single-sided, for original length)
  const halfN = Math.floor(N0 / 2) + 1;
  const amp = new Array(halfN);
  for (let i = 0; i < halfN; i++) {
    amp[i] = 2 * Math.sqrt(real[i] ** 2 + imag[i] ** 2) / N0;
  }
  return amp;
}

/* ========== Harmonic Ratio ========== */
function calcHarmonicRatios(cyc_ap, cyc_ml, cyc_vt) {
  const amp_ap = fft(cyc_ap);
  const amp_ml = fft(cyc_ml);
  const amp_vt = fft(cyc_vt);

  // Harmonics 1-20 (1-indexed)
  let oddAP = 0, evenAP = 0;
  let oddML = 0, evenML = 0;
  let oddVT = 0, evenVT = 0;

  for (let h = 1; h <= 20; h++) {
    if (h >= amp_ap.length) break;
    if (h % 2 === 1) { // odd
      oddAP += amp_ap[h];
      oddML += amp_ml[h];
      oddVT += amp_vt[h];
    } else { // even
      evenAP += amp_ap[h];
      evenML += amp_ml[h];
      evenVT += amp_vt[h];
    }
  }

  return {
    ap: evenAP > 0 ? oddAP / evenAP : 0,  // AP: odd/even
    ml: oddML > 0 ? evenML / oddML : 0,    // ML: even/odd
    vt: evenVT > 0 ? oddVT / evenVT : 0    // VT: odd/even
  };
}

/* ========== Statistics ========== */
function mean(arr) { return arr.reduce((s, v) => s + v, 0) / arr.length; }

/* ========== Linear interpolation for resampling ========== */
function linearInterp(xOld, yOld, xNew) {
  const result = new Array(xNew.length);
  let j = 0;
  for (let i = 0; i < xNew.length; i++) {
    while (j < xOld.length - 2 && xOld[j + 1] < xNew[i]) j++;
    const t = (xNew[i] - xOld[j]) / (xOld[j + 1] - xOld[j]);
    result[i] = yOld[j] + t * (yOld[j + 1] - yOld[j]);
  }
  return result;
}
function sd(arr) {
  const m = mean(arr);
  return Math.sqrt(arr.reduce((s, v) => s + (v - m) ** 2, 0) / (arr.length - 1));
}
function rms(arr) { return Math.sqrt(arr.reduce((s, v) => s + v * v, 0) / arr.length); }

/* ========== Display metrics ========== */
function showMetrics(m) {
  const items = [
    ['歩行速度',        m.walkSpeed.toFixed(3),  'm/s'],
    ['歩行時間',        m.boutDuration.toFixed(2), '秒'],
    ['ストライド時間',  m.strideMean.toFixed(3),  '秒'],
    ['ストライド変動(SD)', m.strideSD.toFixed(4), '秒'],
    ['ストライド変動係数', m.strideCV.toFixed(2),  '%'],
    ['ケイデンス',      m.cadence.toFixed(1),     'steps/min'],
    ['HR 前後(AP)',      m.hrAP.toFixed(3),       ''],
    ['HR 左右(ML)',      m.hrML.toFixed(3),       ''],
    ['HR 上下(VT)',      m.hrVT.toFixed(3),       ''],
    ['RMS 前後(AP)',     m.rmsAP.toFixed(4),      'G'],
    ['RMS 左右(ML)',     m.rmsML.toFixed(4),      'G'],
    ['RMS AP / 速度',   m.rmsAP_norm.toFixed(4),  'G·s/m'],
    ['RMS ML / 速度',   m.rmsML_norm.toFixed(4),  'G·s/m'],
    ['体幹不安定性',     m.trunkInstab.toFixed(4), 'G·s/m'],
    ['有効ストライド数',  m.nStrides,              ''],
  ];

  const el = $('metrics');
  el.innerHTML = items.map(([k, v, u]) =>
    `<div class="metric"><div class="k">${k}</div><div class="v">${v}<span class="u">${u}</span></div></div>`
  ).join('');
  el.style.display = 'grid';
}

/* ========== Plots ========== */
function plotResults(d) {
  const {
    time, filt_ap, filt_ml, filt_vt,
    bout, msd_sum, threshold,
    bout_ap, bout_ml, bout_vt,
    filt_bout_ap, peaks, exceptStride,
    sampleRate,
    strideTimes, HR_ap_arr, HR_ml_arr, HR_vt_arr
  } = d;

  const plotCfg = { responsive: true, displayModeBar: false };
  const layout = (title, xlab, ylab) => ({
    title: { text: title, font: { size: 14 } },
    xaxis: { title: xlab },
    yaxis: { title: ylab },
    margin: { t: 40, b: 40, l: 50, r: 20 },
    legend: { orientation: 'h', y: -0.2 },
    font: { family: 'system-ui, sans-serif' }
  });

  /* Plot 1: Filtered full signal with bout highlight */
  const timeArr = time.map((_, i) => i / sampleRate);
  Plotly.newPlot('plotSignal', [
    { y: filt_ap, x: timeArr, name: 'AP', type: 'scatter', mode: 'lines', line: { width: 1 } },
    { y: filt_ml, x: timeArr, name: 'ML', type: 'scatter', mode: 'lines', line: { width: 1 } },
    { y: filt_vt, x: timeArr, name: 'VT', type: 'scatter', mode: 'lines', line: { width: 1 } },
  ], {
    ...layout('フィルタ後加速度信号', '時間 (秒)', '加速度 (G)'),
    shapes: [{
      type: 'rect',
      x0: bout.start / sampleRate, x1: bout.end / sampleRate,
      y0: 0, y1: 1, yref: 'paper',
      fillcolor: 'rgba(46, 134, 193, 0.15)', line: { width: 0 }
    }],
    annotations: [{
      x: (bout.start + bout.end) / 2 / sampleRate, y: 1, yref: 'paper',
      text: '歩行区間', showarrow: false, font: { size: 12, color: '#2E86C1' }
    }]
  }, plotCfg);

  /* Plot 2: Bout detection (moving SD sum) */
  Plotly.newPlot('plotBout', [
    { y: msd_sum, x: msd_sum.map((_, i) => i / sampleRate), name: '移動SD合計', type: 'scatter', mode: 'lines', line: { width: 1, color: '#555' } },
    { y: msd_sum.map(() => threshold), x: msd_sum.map((_, i) => i / sampleRate), name: `閾値 (${threshold})`, type: 'scatter', mode: 'lines', line: { dash: 'dash', color: 'red', width: 1 } }
  ], layout('歩行区間検出（移動SD）', '時間 (秒)', '移動SD合計'), plotCfg);

  /* Plot 3: Stride detection with peaks */
  const boutTime = bout_ap.map((_, i) => i / sampleRate);
  const peakX = peaks.map(p => p / sampleRate);
  const peakY = peaks.map(p => filt_bout_ap[p]);

  // Mark valid vs excluded peaks
  const nStrides = Math.floor(peaks.length / 2);
  const validStart = exceptStride * 2 - 1;
  const validEnd = (nStrides - exceptStride) * 2;

  Plotly.newPlot('plotStride', [
    { y: bout_ap, x: boutTime, name: 'AP (歩行区間)', type: 'scatter', mode: 'lines', line: { width: 1, color: 'navy' } },
    { y: filt_bout_ap, x: boutTime, name: 'AP (2Hz LPF)', type: 'scatter', mode: 'lines', line: { width: 1.5, color: '#e67e22' } },
    { y: peakY, x: peakX, name: 'ピーク', type: 'scatter', mode: 'markers', marker: { size: 8, color: 'maroon' } },
  ], layout('ストライド検出', '時間 (秒)', '加速度 (G)'), plotCfg);

  /* Plot 4: Per-stride metrics bar chart */
  const strideLabels = strideTimes.map((_, i) => `S${i + 1}`);
  Plotly.newPlot('plotHR', [
    { y: HR_ap_arr, x: strideLabels, name: 'HR AP', type: 'bar' },
    { y: HR_ml_arr, x: strideLabels, name: 'HR ML', type: 'bar' },
    { y: HR_vt_arr, x: strideLabels, name: 'HR VT', type: 'bar' },
  ], {
    ...layout('ストライド別 Harmonic Ratio', 'ストライド', 'HR'),
    barmode: 'group'
  }, plotCfg);

  /* Plot 5: FFT of last valid stride (example) */
  const lastValidP = Math.min((nStrides - exceptStride - 1), nStrides - 1);
  if (lastValidP >= exceptStride && peaks[2 * lastValidP - 1] !== undefined && peaks[2 * lastValidP + 1] !== undefined) {
    const i0 = peaks[2 * lastValidP - 1];
    const i1 = peaks[2 * lastValidP + 1];
    const cyc = bout_ap.slice(i0, i1 + 1);
    const amp = fft(cyc);
    const hz = amp.map((_, i) => i * (sampleRate / cyc.length));

    const cyc_ml_ex = bout_ml.slice(i0, i1 + 1);
    const cyc_vt_ex = bout_vt.slice(i0, i1 + 1);
    const amp_ml_ex = fft(cyc_ml_ex);
    const amp_vt_ex = fft(cyc_vt_ex);

    Plotly.newPlot('plotFFT', [
      { y: amp.slice(0, 21), x: hz.slice(0, 21), name: 'AP', type: 'bar', marker: { color: '#e74c3c' } },
      { y: amp_ml_ex.slice(0, 21), x: hz.slice(0, 21), name: 'ML', type: 'bar', marker: { color: '#27ae60' } },
      { y: amp_vt_ex.slice(0, 21), x: hz.slice(0, 21), name: 'VT', type: 'bar', marker: { color: '#2980b9' } },
    ], {
      ...layout('FFTスペクトル（代表ストライド）', '周波数 (Hz)', '振幅'),
      barmode: 'group'
    }, plotCfg);
  }
}

/* ========== Helpers ========== */
function showMsg(text, cls) {
  const el = $('msg');
  el.textContent = text;
  el.className = cls;
  el.style.display = 'block';
}

/* ========== CSV export ========== */
$('dlCsv').addEventListener('click', () => {
  if (!rawData) return;
  const d = rawData;
  let csv = 'metric,value,unit\n';
  csv += `walk_speed,${d.walkSpeed.toFixed(4)},m/s\n`;
  csv += `bout_duration,${d.boutDuration.toFixed(3)},s\n`;
  csv += `stride_mean,${d.strideMean.toFixed(4)},s\n`;
  csv += `stride_sd,${d.strideSD.toFixed(5)},s\n`;
  csv += `stride_cv,${d.strideCV.toFixed(3)},%\n`;
  csv += `cadence,${d.cadence.toFixed(2)},steps/min\n`;
  csv += `HR_AP,${d.hrAP.toFixed(4)},\n`;
  csv += `HR_ML,${d.hrML.toFixed(4)},\n`;
  csv += `HR_VT,${d.hrVT.toFixed(4)},\n`;
  csv += `RMS_AP,${d.rmsAP.toFixed(5)},G\n`;
  csv += `RMS_ML,${d.rmsML.toFixed(5)},G\n`;
  csv += `RMS_AP_norm,${d.rmsAP_norm.toFixed(5)},G·s/m\n`;
  csv += `RMS_ML_norm,${d.rmsML_norm.toFixed(5)},G·s/m\n`;
  csv += `trunk_instability,${d.trunkInstab.toFixed(5)},G·s/m\n`;
  csv += `n_strides,${d.nStrides},\n`;

  csv += '\nstride,stride_time,HR_AP,HR_ML,HR_VT,RMS_AP,RMS_ML\n';
  for (let i = 0; i < d.strideTimes.length; i++) {
    csv += `${i + 1},${d.strideTimes[i].toFixed(4)},${d.HR_ap_arr[i].toFixed(4)},${d.HR_ml_arr[i].toFixed(4)},${d.HR_vt_arr[i].toFixed(4)},${d.rms_ap_arr[i].toFixed(5)},${d.rms_ml_arr[i].toFixed(5)}\n`;
  }

  download(csv, 'gait_result.csv', 'text/csv');
});

/* ========== PNG export ========== */
$('dlPng').addEventListener('click', async () => {
  const plots = ['plotSignal', 'plotBout', 'plotStride', 'plotHR', 'plotFFT'];
  for (const id of plots) {
    const el = document.getElementById(id);
    if (el && el.data) {
      await Plotly.downloadImage(el, { format: 'png', width: 1200, height: 500, filename: id });
    }
  }
});

function download(content, filename, type) {
  const blob = new Blob([content], { type });
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = filename;
  a.click();
  URL.revokeObjectURL(a.href);
}

/* ========== Populate column selects on file change ========== */
$('file').addEventListener('change', async function () {
  const file = this.files[0];
  if (!file) return;
  const text = await file.text();
  const first = Papa.parse(text.trim(), { header: true, preview: 2 });
  const fields = first.meta.fields || [];
  const selects = ['colAP', 'colML', 'colVT'];
  const defaults = [/gfx|ap|anteroposterior/i, /gfy|ml|mediolateral/i, /gfz|vt|vertical/i];
  selects.forEach((id, si) => {
    const sel = $(id);
    sel.innerHTML = '<option value="">自動検出</option>';
    fields.forEach(f => {
      const opt = document.createElement('option');
      opt.value = f;
      opt.textContent = f;
      if (defaults[si].test(f)) opt.selected = true;
      sel.appendChild(opt);
    });
  });
});
