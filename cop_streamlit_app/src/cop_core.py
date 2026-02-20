"""
COP (Center of Pressure) analysis utilities.

Expected CSV columns (case-insensitive, surrounding spaces ignored):
- time : seconds
- COPx : medio-lateral (ML)
- COPy : antero-posterior (AP)

The app trims the first/last `trim_seconds`, centers COP, low-pass filters,
converts to cm, then computes:
- 95% confidence ellipse area (cm^2)
- convex hull (outer) area (cm^2)
- band power ratios (0–3 Hz normalized to 100%)
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Tuple, Optional

import numpy as np
import pandas as pd
from scipy.signal import butter, filtfilt, welch
from scipy.spatial import ConvexHull
from scipy.stats import chi2


BANDS_DEFAULT: Dict[str, Tuple[float, float]] = {
    "LF": (0.0, 0.3),
    "MF": (0.3, 1.0),
    "HF": (1.0, 3.0),
}


@dataclass(frozen=True)
class CopConfig:
    fs: float = 100.0          # sampling frequency (Hz)
    lowpass_hz: float = 10.0   # low-pass cutoff (Hz)
    trim_seconds: float = 5.0  # trim first/last seconds
    unit_to_cm: float = 0.1    # mm->cm: 0.1, already cm: 1.0
    bands: Dict[str, Tuple[float, float]] = None  # defaults to BANDS_DEFAULT

    def __post_init__(self):
        # dataclass frozen workaround: validate via object.__setattr__ not allowed; keep light validation.
        pass


@dataclass
class CopResult:
    # time series after trimming/centering/filtering (cm)
    time: np.ndarray
    copx_cm: np.ndarray
    copy_cm: np.ndarray
    copx_f: np.ndarray  # filtered (original units, before cm conversion)
    copy_f: np.ndarray

    # summary metrics
    ellipse_area_cm2: float
    hull_area_cm2: float

    # band ratios (%), per axis and combined
    band_ratios_ml: Dict[str, float]
    band_ratios_ap: Dict[str, float]
    band_ratios_total: Dict[str, float]

    # confidence ellipse geometry
    ellipse_xy: np.ndarray  # shape (2, N) centered at (0,0)


def _norm_col(s: str) -> str:
    return str(s).strip().lower().replace(" ", "")


def lowpass_filter(x: np.ndarray, fs: float, cutoff: float, order: int = 4) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    # filtfilt requires len(x) > padlen; use a conservative threshold
    if len(x) < (order * 3 + 15):
        return x
    nyq = 0.5 * fs
    if cutoff <= 0 or cutoff >= nyq:
        return x
    b, a = butter(order, cutoff / nyq, btype="low")
    return filtfilt(b, a, x)


def _band_ratios_1d(x: np.ndarray, fs: float, bands: Dict[str, Tuple[float, float]]) -> Dict[str, float]:
    x = np.asarray(x, dtype=float)
    nperseg = int(min(1024, max(16, len(x))))
    f, pxx = welch(x, fs=fs, nperseg=nperseg)

    # Normalize by total power in 0–3 Hz
    total_mask = (f >= 0.0) & (f <= 3.0)
    total_power = float(np.trapz(pxx[total_mask], f[total_mask])) if np.any(total_mask) else 0.0

    out: Dict[str, float] = {}
    for name, (f1, f2) in bands.items():
        mask = (f >= f1) & (f <= f2)
        if np.any(mask) and total_power > 0:
            band_power = float(np.trapz(pxx[mask], f[mask]))
            out[name] = band_power / total_power * 100.0
        else:
            out[name] = float("nan")
    return out


def _band_ratios_total(copx_f: np.ndarray, copy_f: np.ndarray, fs: float, bands: Dict[str, Tuple[float, float]]) -> Dict[str, float]:
    nperseg = int(min(1024, max(16, len(copx_f))))
    f, pxx = welch(copx_f, fs=fs, nperseg=nperseg)
    _, pyy = welch(copy_f, fs=fs, nperseg=nperseg)
    p_total = pxx + pyy

    total_mask = (f >= 0.0) & (f <= 3.0)
    total_power = float(np.trapz(p_total[total_mask], f[total_mask])) if np.any(total_mask) else 0.0

    out: Dict[str, float] = {}
    for name, (f1, f2) in bands.items():
        mask = (f >= f1) & (f <= f2)
        if np.any(mask) and total_power > 0:
            band_power = float(np.trapz(p_total[mask], f[mask]))
            out[name] = band_power / total_power * 100.0
        else:
            out[name] = float("nan")
    return out


def _confidence_ellipse_xy(copx_cm: np.ndarray, copy_cm: np.ndarray, conf: float = 0.95, n: int = 400) -> np.ndarray:
    cov = np.cov(copx_cm, copy_cm)
    eigvals, eigvecs = np.linalg.eigh(cov)
    order = np.argsort(eigvals)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]

    scale = float(np.sqrt(chi2.ppf(conf, df=2)))
    a = scale * np.sqrt(max(eigvals[0], 0.0))
    b = scale * np.sqrt(max(eigvals[1], 0.0))
    angle = float(np.arctan2(eigvecs[1, 0], eigvecs[0, 0]))

    t = np.linspace(0, 2 * np.pi, n)
    ellipse = np.array([a * np.cos(t), b * np.sin(t)])
    R = np.array([[np.cos(angle), -np.sin(angle)],
                  [np.sin(angle),  np.cos(angle)]])
    return R @ ellipse


def analyze_cop_csv(csv_path: str | Path, cfg: Optional[CopConfig] = None) -> CopResult:
    """
    Load a COP CSV and compute metrics + ellipse geometry.

    Raises ValueError if required columns are missing.
    """
    if cfg is None:
        cfg = CopConfig()
    bands = cfg.bands or BANDS_DEFAULT

    csv_path = Path(csv_path)
    df = pd.read_csv(csv_path)
    df.columns = [str(c).strip() for c in df.columns]
    colmap = {_norm_col(c): c for c in df.columns}

    required = {"time", "copx", "copy"}
    if not required.issubset(set(colmap.keys())):
        raise ValueError(f"CSV must contain columns: time, COPx, COPy (found: {list(df.columns)})")

    time = pd.to_numeric(df[colmap["time"]], errors="coerce").to_numpy(dtype=float)
    copx = pd.to_numeric(df[colmap["copx"]], errors="coerce").to_numpy(dtype=float)
    copy = pd.to_numeric(df[colmap["copy"]], errors="coerce").to_numpy(dtype=float)

    # drop NaN rows
    mask = ~np.isnan(time) & ~np.isnan(copx) & ~np.isnan(copy)
    time, copx, copy = time[mask], copx[mask], copy[mask]

    # sort by time (just in case)
    if len(time) > 1:
        order = np.argsort(time)
        time, copx, copy = time[order], copx[order], copy[order]

    # trim first/last seconds
    if len(time) > 0 and cfg.trim_seconds > 0:
        t0 = time[0] + cfg.trim_seconds
        t1 = time[-1] - cfg.trim_seconds
        mask_t = (time >= t0) & (time <= t1)
        if np.any(mask_t):
            time, copx, copy = time[mask_t], copx[mask_t], copy[mask_t]

    # center
    copx = copx - np.mean(copx) if len(copx) else copx
    copy = copy - np.mean(copy) if len(copy) else copy

    # filter (in original units)
    copx_f = lowpass_filter(copx, fs=cfg.fs, cutoff=cfg.lowpass_hz)
    copy_f = lowpass_filter(copy, fs=cfg.fs, cutoff=cfg.lowpass_hz)

    # unit conversion to cm
    copx_cm = copx_f * cfg.unit_to_cm
    copy_cm = copy_f * cfg.unit_to_cm

    # ellipse area (cm^2)
    cov = np.cov(copx_cm, copy_cm)
    eigvals = np.linalg.eigvalsh(cov)  # ascending
    scale = float(np.sqrt(chi2.ppf(0.95, df=2)))
    ellipse_area_cm2 = float(np.pi * (scale**2) * np.sqrt(max(eigvals[0], 0.0) * max(eigvals[1], 0.0)))

    # hull (outer) area (cm^2)
    points_cm = np.column_stack([copx_cm, copy_cm])
    hull_area_cm2 = float("nan")
    if len(points_cm) >= 3:
        hull = ConvexHull(points_cm)
        hull_area_cm2 = float(hull.volume)  # 2D: area

    # band ratios
    band_ratios_ml = _band_ratios_1d(copx_f, cfg.fs, bands=bands)
    band_ratios_ap = _band_ratios_1d(copy_f, cfg.fs, bands=bands)
    band_ratios_total = _band_ratios_total(copx_f, copy_f, cfg.fs, bands=bands)

    ellipse_xy = _confidence_ellipse_xy(copx_cm, copy_cm, conf=0.95, n=400)

    return CopResult(
        time=time,
        copx_cm=copx_cm,
        copy_cm=copy_cm,
        copx_f=copx_f,
        copy_f=copy_f,
        ellipse_area_cm2=ellipse_area_cm2,
        hull_area_cm2=hull_area_cm2,
        band_ratios_ml=band_ratios_ml,
        band_ratios_ap=band_ratios_ap,
        band_ratios_total=band_ratios_total,
        ellipse_xy=ellipse_xy,
    )
