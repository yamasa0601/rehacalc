
import io
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import streamlit as st

from emg_processing import parse_emg_csv, analyze_signal

st.set_page_config(page_title="EMG Heel-Strike Estimator", layout="wide")

st.title("EMGからHeel Strike推定 → %peak → 歩行周期(0–100%)プロット")
st.caption("歩行のみ想定。50–450Hz bandpass → 整流 → 50ms RMS → バースト立ち上がり=HS（動的閾値 median + k×MAD）。")

with st.sidebar:
    st.header("設定（基本）")
    fs = st.number_input("サンプリング周波数 (Hz)", min_value=100, max_value=5000, value=1000, step=100)

    st.header("前処理")
    hp = st.number_input("HighPass (Hz)", min_value=1.0, max_value=500.0, value=50.0, step=1.0)
    lp = st.number_input("LowPass (Hz)",  min_value=10.0, max_value=1000.0, value=450.0, step=10.0)
    rms_ms = st.number_input("RMS窓 (ms)", min_value=1.0, max_value=500.0, value=50.0, step=1.0)

    st.header("HS推定（バースト）")
    min_burst = st.number_input("最小バースト持続 (ms)", min_value=5.0, max_value=500.0, value=30.0, step=5.0)
    min_gap = st.number_input("最小HS間隔 (ms)（二重検出抑制）", min_value=50.0, max_value=2000.0, value=300.0, step=50.0)
    ic_offset = st.number_input("HSオフセット (ms)（onset→HS補正）", min_value=-200.0, max_value=200.0, value=0.0, step=5.0)

    st.header("歩数（任意）")
    total_steps = st.number_input("動画の歩数（分かるなら）", min_value=0, max_value=500, value=0, step=1)
    st.caption("入力すると、推定HS数が近くなるように自動調整します（左右で12/13のように分配）。")

    st.header("出力")
    n_points = st.number_input("周期正規化の点数", min_value=101, max_value=2001, value=501, step=50)
    save_fig = st.checkbox("図を保存（PNG）", value=True)
    dpi = st.number_input("保存DPI", min_value=72, max_value=600, value=300, step=10)

st.subheader("CSVアップロード（2ch想定：左/右など）")
col_u1, col_u2 = st.columns(2)

with col_u1:
    f1 = st.file_uploader("CSV 1（例：Lt）", type=["csv"], accept_multiple_files=False)
with col_u2:
    f2 = st.file_uploader("CSV 2（例：Rt）", type=["csv"], accept_multiple_files=False)

run = st.button("解析してプロット", type="primary", use_container_width=True)

def _process_one(uploaded_file, expected_hs=None):
    b = uploaded_file.getvalue()
    sig = parse_emg_csv(b, filename=uploaded_file.name, fs=float(fs))
    res = analyze_signal(
        sig,
        highpass_hz=float(hp),
        lowpass_hz=float(lp),
        rms_ms=float(rms_ms),
        expected_hs=expected_hs,
        min_burst_ms=float(min_burst),
        min_gap_ms=float(min_gap),
        ic_offset_ms=float(ic_offset),
        n_points=int(n_points),
    )
    return sig, res

def _plot_debug(sig, res):
    t = sig.t
    hs = res.hs_idx
    fig = plt.figure(figsize=(10, 4))
    plt.plot(t, res.env_pctpeak, label="%peak envelope")
    plt.axhline((res.threshold / np.nanmax(res.env) * 100.0) if np.nanmax(res.env) > 0 else np.nan, linestyle="--", label="threshold (scaled)")
    for idx in hs:
        plt.axvline(t[int(idx)], alpha=0.25)
    plt.title(f"{res.name}: %peak envelope + estimated HS (n={len(hs)})")
    plt.xlabel("time (s)")
    plt.ylabel("%RMS (%peak)")
    plt.grid(True, alpha=0.25)
    plt.tight_layout()
    return fig

def _plot_cycles(res):
    grid = res.meta["grid"]
    cycles = res.meta["cycles"]
    fig = plt.figure(figsize=(10, 5))
    if cycles.shape[0] > 0:
        # individual
        for i in range(cycles.shape[0]):
            plt.plot(grid, cycles[i], color="lightgray", alpha=0.35, linewidth=0.7)
        mean = np.nanmean(cycles, axis=0)
        sd = np.nanstd(cycles, axis=0)
        plt.plot(grid, mean, linewidth=2.2, label="mean")
        plt.fill_between(grid, mean - sd, mean + sd, alpha=0.2, label="±1SD")
    plt.title(f"{res.name}: gait cycle normalized (0–100%)  %peak")
    plt.xlabel("gait cycle (%)")
    plt.ylabel("%RMS (%peak)")
    plt.xlim(0, 100)
    plt.grid(True, alpha=0.25)
    plt.tight_layout()
    return fig

def _save_figure(fig, out_dir: Path, stem: str, suffix: str, dpi: int):
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"{stem}_{suffix}.png"
    fig.savefig(path, dpi=dpi)
    return path

if run:
    if f1 is None or f2 is None:
        st.error("CSVを2つアップロードしてください。")
        st.stop()

    exp1 = exp2 = None
    if int(total_steps) > 0:
        exp1 = int(total_steps) // 2
        exp2 = int(total_steps) - exp1

    sig1, res1 = _process_one(f1, expected_hs=exp1)
    sig2, res2 = _process_one(f2, expected_hs=exp2)

    out_dir = Path("outputs")

    # Summary table
    st.subheader("推定結果サマリ")
    sum_df = pd.DataFrame([
        {"signal": res1.name, **{k: v for k, v in res1.stats.items()}, "k_mad": res1.k_mad},
        {"signal": res2.name, **{k: v for k, v in res2.stats.items()}, "k_mad": res2.k_mad},
    ])
    st.dataframe(sum_df, use_container_width=True)

    # HS list download
    st.subheader("推定HS（サンプル番号と時刻）")
    def hs_table(sig, res):
        hs_no = (res.hs_idx.astype(int) + 1)
        hs_t = sig.t[res.hs_idx.astype(int)]
        return pd.DataFrame({"hs_index": np.arange(1, len(hs_no)+1), "No": hs_no, "time_s": hs_t})

    hs1 = hs_table(sig1, res1)
    hs2 = hs_table(sig2, res2)
    col_h1, col_h2 = st.columns(2)
    with col_h1:
        st.write(f"**{res1.name}**")
        st.dataframe(hs1, use_container_width=True, height=260)
        st.download_button("HS CSVをダウンロード", data=hs1.to_csv(index=False).encode("utf-8-sig"),
                           file_name=f"{res1.name}_HS.csv", mime="text/csv")
    with col_h2:
        st.write(f"**{res2.name}**")
        st.dataframe(hs2, use_container_width=True, height=260)
        st.download_button("HS CSVをダウンロード", data=hs2.to_csv(index=False).encode("utf-8-sig"),
                           file_name=f"{res2.name}_HS.csv", mime="text/csv")

    st.subheader("プロット")
    col_p1, col_p2 = st.columns(2)

    with col_p1:
        st.write(f"**{res1.name}**")
        fig_debug_1 = _plot_debug(sig1, res1)
        st.pyplot(fig_debug_1, use_container_width=True)
        fig_cycle_1 = _plot_cycles(res1)
        st.pyplot(fig_cycle_1, use_container_width=True)

        if save_fig:
            p1a = _save_figure(fig_debug_1, out_dir, res1.name, "debug_envelope_hs", int(dpi))
            p1b = _save_figure(fig_cycle_1, out_dir, res1.name, "cycle_mean_sd", int(dpi))
            st.success(f"Saved: {p1a} / {p1b}")

    with col_p2:
        st.write(f"**{res2.name}**")
        fig_debug_2 = _plot_debug(sig2, res2)
        st.pyplot(fig_debug_2, use_container_width=True)
        fig_cycle_2 = _plot_cycles(res2)
        st.pyplot(fig_cycle_2, use_container_width=True)

        if save_fig:
            p2a = _save_figure(fig_debug_2, out_dir, res2.name, "debug_envelope_hs", int(dpi))
            p2b = _save_figure(fig_cycle_2, out_dir, res2.name, "cycle_mean_sd", int(dpi))
            st.success(f"Saved: {p2a} / {p2b}")

    st.info("※HS推定は“筋活動バーストの立ち上がり”を利用するため、筋の種類・歩容・貼付位置によってはHSとズレます。必要ならIC_OFFSET(ms)で補正してください。")
