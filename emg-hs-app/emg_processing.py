
from __future__ import annotations

import io
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Tuple, Dict, Any, List, Union

import numpy as np
import pandas as pd
from scipy.signal import butter, sosfiltfilt


@dataclass
class Signal:
    name: str
    fs: float
    t: np.ndarray          # seconds
    x: np.ndarray          # raw signal (mV or a.u.)
    meta: Dict[str, Any]


@dataclass
class HSResult:
    name: str
    fs: float
    hs_idx: np.ndarray     # sample indices (0-based)
    k_mad: float
    threshold: float
    env: np.ndarray        # envelope (RMS of rectified)
    env_pctpeak: np.ndarray
    stats: Dict[str, float]
    meta: Dict[str, Any]


def _mad(x: np.ndarray) -> float:
    med = np.nanmedian(x)
    return float(np.nanmedian(np.abs(x - med)))


def _rising_edges(mask: np.ndarray):
    """Return list of (start, end) for True segments [start, end)."""
    mask = mask.astype(bool)
    d = np.diff(mask.astype(int))
    starts = np.where(d == 1)[0] + 1
    ends   = np.where(d == -1)[0] + 1
    if mask[0]:
        starts = np.r_[0, starts]
    if mask[-1]:
        ends = np.r_[ends, len(mask)]
    return list(zip(starts, ends))


def bandpass_filter(x: np.ndarray, fs: float, lowcut: float, highcut: float, order: int = 4) -> np.ndarray:
    sos = butter(order, [lowcut, highcut], btype="bandpass", fs=fs, output="sos")
    return sosfiltfilt(sos, x)


def moving_rms(x: np.ndarray, win: int) -> np.ndarray:
    if win <= 1:
        return np.abs(x)
    w = np.ones(win, dtype=float) / win
    return np.sqrt(np.convolve(x * x, w, mode="same"))


def parse_emg_csv(
    file_bytes: bytes,
    filename: str,
    fs: float = 1000.0,
) -> Signal:
    """
    Parse either:
      A) "sensor style" CSV with meta header then a table starting with 'No,...' (your example), or
      B) a normal table CSV with columns like: time, ch1 (and maybe others).

    Returns a Signal with time in seconds and x as the main channel signal.
    """
    # Try decode (cp932 first, then utf-8-sig)
    for enc in ("cp932", "utf-8-sig", "utf-8"):
        try:
            text = file_bytes.decode(enc)
            used_enc = enc
            break
        except Exception:
            continue
    else:
        text = file_bytes.decode("utf-8", errors="replace")
        used_enc = "utf-8 (errors=replace)"

    lines = text.splitlines()

    # Detect sensor style: find first line starting with "No,"
    start = None
    for i, line in enumerate(lines):
        if line.startswith("No,") or line.startswith("No\t"):
            start = i
            break

    meta: Dict[str, Any] = {"encoding": used_enc, "source_filename": filename}

    if start is not None:
        # Sensor style
        # Try to extract a human name
        sensor_name = Path(filename).stem
        for l in lines[:80]:
            if l.startswith("センサ名"):
                parts = [p.strip() for p in l.split(",")]
                if len(parts) >= 2 and parts[1]:
                    sensor_name = parts[1]
                    break
        meta["sensor_name"] = sensor_name

        # If fs appears, read it, else trust passed fs
        for l in lines[:120]:
            if "サンプリング周波数" in l:
                parts = [p.strip() for p in l.split(",") if p.strip()]
                for p in parts:
                    if re.fullmatch(r"\d+(\.\d+)?", p):
                        fs = float(p)
                        meta["fs_from_file"] = fs
                        break

        df = pd.read_csv(io.StringIO("\n".join(lines[start:])))
        df.columns = [c.strip() for c in df.columns]
        if "No" not in df.columns:
            raise ValueError("Could not find 'No' column after header.")
        # voltage column
        vcol = None
        for c in df.columns:
            if "電圧" in c or "mV" in c or "voltage" in c.lower():
                vcol = c
                break
        if vcol is None:
            # fallback: second column
            vcol = df.columns[1]

        no = df["No"].astype(int).to_numpy()
        x = pd.to_numeric(df[vcol], errors="coerce").to_numpy(dtype=float)
        t = (no - 1) / fs
        name = sensor_name

        return Signal(name=name, fs=fs, t=t, x=x, meta={**meta, "format": "sensor_style", "value_col": vcol})

    # Otherwise normal CSV table
    df = pd.read_csv(io.BytesIO(file_bytes))
    df.columns = [str(c).strip() for c in df.columns]

    name = Path(filename).stem
    time_col = None
    for c in df.columns:
        if c.lower() in ("time", "time_s", "t", "sec", "seconds"):
            time_col = c
            break

    if time_col is not None:
        t = pd.to_numeric(df[time_col], errors="coerce").to_numpy(dtype=float)
        # Choose first non-time column as signal
        sig_col = [c for c in df.columns if c != time_col][0]
        x = pd.to_numeric(df[sig_col], errors="coerce").to_numpy(dtype=float)
        return Signal(name=name, fs=fs, t=t, x=x, meta={**meta, "format": "table", "time_col": time_col, "value_col": sig_col})

    # If no time, assume first column is sample index
    idx_col = df.columns[0]
    x_col = df.columns[1] if len(df.columns) > 1 else df.columns[0]
    idx = pd.to_numeric(df[idx_col], errors="coerce").to_numpy(dtype=float)
    x = pd.to_numeric(df[x_col], errors="coerce").to_numpy(dtype=float)
    t = (idx - idx[0]) / fs
    return Signal(name=name, fs=fs, t=t, x=x, meta={**meta, "format": "table_no_time", "index_col": idx_col, "value_col": x_col})


def make_envelope(
    x: np.ndarray,
    fs: float,
    highpass_hz: float = 50.0,
    lowpass_hz: float = 450.0,
    rms_ms: float = 50.0,
    order: int = 4,
) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    x = x - np.nanmean(x)
    bp = bandpass_filter(x, fs, highpass_hz, lowpass_hz, order=order)
    rect = np.abs(bp)
    win = int(round(fs * rms_ms / 1000.0))
    env = moving_rms(rect, win)
    return env


def percent_peak(env: np.ndarray, hs_idx: Optional[np.ndarray] = None) -> np.ndarray:
    env = np.asarray(env, dtype=float)
    if hs_idx is not None and len(hs_idx) >= 2:
        s = int(hs_idx[0]); e = int(hs_idx[-1])
        ref = np.nanmax(env[s:e])
    else:
        ref = np.nanmax(env)
    if not np.isfinite(ref) or ref <= 0:
        return np.full_like(env, np.nan)
    return (env / ref) * 100.0


def detect_hs_from_bursts(
    env: np.ndarray,
    fs: float,
    expected_count: Optional[int] = None,
    k_candidates: Optional[np.ndarray] = None,
    min_burst_ms: float = 30.0,
    min_gap_ms: float = 300.0,
    ic_offset_ms: float = 0.0,
    plausible_step_s: Tuple[float, float] = (0.35, 1.6),
) -> Tuple[float, float, np.ndarray, Dict[str, float]]:
    """
    Return (k_mad, threshold, hs_idx, stats).
    hs_idx are 0-based indices (burst onsets merged by min_gap, plus offset).
    """
    env = np.asarray(env, dtype=float)

    if k_candidates is None:
        k_candidates = np.round(np.arange(1.2, 4.1, 0.1), 2)

    med = np.nanmedian(env)
    m = _mad(env)
    min_burst = int(round(fs * min_burst_ms / 1000.0))
    min_gap = int(round(fs * min_gap_ms / 1000.0))
    offset = int(round(fs * ic_offset_ms / 1000.0))

    best_score = None
    best_pack = None

    for k in k_candidates:
        thr = med + float(k) * m
        bursts = _rising_edges(env > thr)
        bursts = [(s, e) for (s, e) in bursts if (e - s) >= min_burst]

        hs = np.array([s for (s, e) in bursts], dtype=int)
        if len(hs) == 0:
            continue
        hs = np.sort(hs)

        # merge close detections
        merged = [int(hs[0])]
        for idx in hs[1:]:
            if int(idx) - merged[-1] >= min_gap:
                merged.append(int(idx))
        hs = np.array(merged, dtype=int) + offset
        hs = hs[(hs >= 0) & (hs < len(env))]

        if len(hs) < 2:
            continue

        # scoring
        intervals = np.diff(hs) / fs
        med_int = float(np.median(intervals))
        cv_int = float(np.std(intervals) / (np.mean(intervals) + 1e-12))

        # penalty for implausible intervals
        lo, hi = plausible_step_s
        bad = np.mean((intervals < lo) | (intervals > hi))
        penalty = 5.0 * bad + 2.0 * cv_int + 0.5 * abs(med_int - 0.7)  # 0.7s ~ moderate walking

        if expected_count is not None:
            penalty += 0.8 * abs(len(hs) - expected_count)

        score = penalty

        if best_score is None or score < best_score:
            best_score = score
            best_pack = (float(k), float(thr), hs, {"n_hs": float(len(hs)), "median_step_s": med_int, "cv_step": cv_int, "bad_interval_ratio": float(bad)})

    if best_pack is None:
        raise ValueError("HS detection failed. Try lowering k range or min_gap_ms, or check signal quality.")

    return best_pack


def normalize_cycles(
    y: np.ndarray,
    hs_idx: np.ndarray,
    t: np.ndarray,
    n_points: int = 501,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Return (percent_grid, cycles_matrix) where cycles_matrix shape = (n_cycles, n_points)
    """
    hs_idx = np.asarray(hs_idx, dtype=int)
    y = np.asarray(y, dtype=float)
    t = np.asarray(t, dtype=float)

    grid = np.linspace(0.0, 100.0, n_points)
    cycles: List[np.ndarray] = []

    for i in range(len(hs_idx) - 1):
        s = int(hs_idx[i]); e = int(hs_idx[i+1])
        if e <= s + 1:
            continue
        seg_t = t[s:e]
        seg_y = y[s:e]

        dur = seg_t[-1] - seg_t[0]
        if not np.isfinite(dur) or dur <= 0:
            continue
        pct = (seg_t - seg_t[0]) / dur * 100.0

        # safe sort unique
        idx = np.argsort(pct)
        x = pct[idx]; v = seg_y[idx]
        uniq = np.concatenate(([True], np.diff(x) > 1e-10))
        x = x[uniq]; v = v[uniq]
        if len(x) < 2:
            continue

        interp = np.interp(grid, x, v)
        cycles.append(interp)

    if len(cycles) == 0:
        return grid, np.empty((0, n_points), dtype=float)

    return grid, np.vstack(cycles)


def analyze_signal(
    sig: Signal,
    highpass_hz: float = 50.0,
    lowpass_hz: float = 450.0,
    rms_ms: float = 50.0,
    expected_hs: Optional[int] = None,
    min_burst_ms: float = 30.0,
    min_gap_ms: float = 300.0,
    ic_offset_ms: float = 0.0,
    n_points: int = 501,
) -> HSResult:
    env = make_envelope(sig.x, sig.fs, highpass_hz, lowpass_hz, rms_ms)
    k_mad, thr, hs_idx, stats = detect_hs_from_bursts(
        env, sig.fs,
        expected_count=expected_hs,
        min_burst_ms=min_burst_ms,
        min_gap_ms=min_gap_ms,
        ic_offset_ms=ic_offset_ms,
    )
    env_pct = percent_peak(env, hs_idx=hs_idx)

    grid, cycles = normalize_cycles(env_pct, hs_idx, sig.t, n_points=n_points)

    # Additional stats
    out_stats = dict(stats)
    out_stats["n_cycles"] = float(max(0, len(hs_idx) - 1))
    if cycles.shape[0] > 0:
        out_stats["mean_of_mean_pct"] = float(np.nanmean(np.nanmean(cycles, axis=1)))
    else:
        out_stats["mean_of_mean_pct"] = float("nan")

    return HSResult(
        name=sig.name,
        fs=sig.fs,
        hs_idx=hs_idx,
        k_mad=k_mad,
        threshold=thr,
        env=env,
        env_pctpeak=env_pct,
        stats=out_stats,
        meta={"signal_meta": sig.meta, "grid": grid, "cycles": cycles},
    )
