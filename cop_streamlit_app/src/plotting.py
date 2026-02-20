from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
from typing import Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull

from .cop_core import CopResult


def make_cop_figure(res: CopResult, xlim: Tuple[float, float] = (-6, 6), ylim: Tuple[float, float] = (-6, 6)) -> plt.Figure:
    """COP trajectory with convex hull + 95% confidence ellipse."""
    points = np.column_stack([res.copx_cm, res.copy_cm])
    hull = ConvexHull(points) if len(points) >= 3 else None

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.plot(res.copx_cm, res.copy_cm, linewidth=2, alpha=0.9)

    if hull is not None:
        for simplex in hull.simplices:
            ax.plot(points[simplex, 0], points[simplex, 1], linewidth=2)

    ax.plot(res.ellipse_xy[0], res.ellipse_xy[1], linewidth=3, linestyle=(0, (8, 6)))
    ax.scatter(0, 0, s=120, facecolors="white", edgecolors="black", linewidths=2)

    ax.set_title(f"95% Confidence Ellipse (Area: {res.ellipse_area_cm2:.2f} cm$^2$)")
    ax.set_xlabel("ML (cm)")
    ax.set_ylabel("AP (cm)")
    ax.grid(True, alpha=0.4)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    fig.tight_layout()
    return fig


def make_band_bar_figure(res: CopResult) -> plt.Figure:
    """Band power ratios bar chart (ML vs AP)."""
    labels = list(res.band_ratios_ml.keys())
    ml_vals = [res.band_ratios_ml[k] for k in labels]
    ap_vals = [res.band_ratios_ap[k] for k in labels]

    x = np.arange(len(labels))
    width = 0.35

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(x - width/2, ml_vals, width, label="ML")
    ax.bar(x + width/2, ap_vals, width, label="AP")
    ax.set_ylabel("Power Ratio (%)")
    ax.set_ylim(0, 100)
    ax.set_xticks(x, labels)
    ax.set_title("Band Power Ratios (0–3 Hz = 100%)")
    ax.grid(axis="y", linestyle="--", alpha=0.5)
    ax.legend()
    fig.tight_layout()
    return fig


def make_panel_figure(res: CopResult, xlim: Tuple[float, float] = (-6, 6), ylim: Tuple[float, float] = (-6, 6)) -> plt.Figure:
    """Two-panel figure: COP+ellipse+hull and band bars."""
    points = np.column_stack([res.copx_cm, res.copy_cm])
    hull = ConvexHull(points) if len(points) >= 3 else None

    labels = list(res.band_ratios_ml.keys())
    ml_vals = [res.band_ratios_ml[k] for k in labels]
    ap_vals = [res.band_ratios_ap[k] for k in labels]

    x = np.arange(len(labels))
    width = 0.35

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    ax = axes[0]
    ax.plot(res.copx_cm, res.copy_cm, linewidth=2, alpha=0.9)
    if hull is not None:
        for simplex in hull.simplices:
            ax.plot(points[simplex, 0], points[simplex, 1], linewidth=2)
    ax.plot(res.ellipse_xy[0], res.ellipse_xy[1], linewidth=3, linestyle=(0, (8, 6)))
    ax.scatter(0, 0, s=120, facecolors="white", edgecolors="black", linewidths=2)
    ax.set_title(f"95% Confidence Ellipse (Area: {res.ellipse_area_cm2:.2f} cm$^2$)")
    ax.set_xlabel("ML (cm)")
    ax.set_ylabel("AP (cm)")
    ax.grid(True, alpha=0.4)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)

    ax = axes[1]
    ax.bar(x - width/2, ml_vals, width, label="ML")
    ax.bar(x + width/2, ap_vals, width, label="AP")
    ax.set_ylabel("Power Ratio (%)")
    ax.set_ylim(0, 100)
    ax.set_xticks(x, labels)
    ax.set_title("Band Power Ratios (0–3 Hz = 100%)")
    ax.grid(axis="y", linestyle="--", alpha=0.5)
    ax.legend()

    fig.tight_layout()
    return fig


def result_tables(res: CopResult) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (metrics_df, bands_df)."""
    metrics = pd.DataFrame(
        {
            "metric": ["ellipse_area_cm2", "hull_area_cm2"],
            "value": [res.ellipse_area_cm2, res.hull_area_cm2],
        }
    )

    bands = pd.DataFrame(
        {
            "band": list(res.band_ratios_ml.keys()),
            "ML_%": [res.band_ratios_ml[k] for k in res.band_ratios_ml.keys()],
            "AP_%": [res.band_ratios_ap[k] for k in res.band_ratios_ap.keys()],
            "TOTAL_%": [res.band_ratios_total[k] for k in res.band_ratios_total.keys()],
        }
    )
    return metrics, bands
