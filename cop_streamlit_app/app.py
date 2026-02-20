import io
from datetime import datetime
from pathlib import Path

import pandas as pd
import streamlit as st

from src.cop_core import CopConfig, analyze_cop_csv
from src.plotting import make_panel_figure, result_tables

APP_DIR = Path(__file__).resolve().parent
DATA_DIR = APP_DIR / "data"
UPLOAD_DIR = DATA_DIR / "uploads"
OUTPUT_DIR = APP_DIR / "outputs"

UPLOAD_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


st.set_page_config(page_title="COP Analysis", layout="wide")
st.title("COP Analysis (CSVを読み込むだけで可視化)")

with st.sidebar:
    st.header("設定")
    fs = st.number_input("サンプリング周波数 FS (Hz)", min_value=1.0, value=100.0, step=1.0)
    lowpass = st.number_input("ローパス cutoff (Hz)", min_value=0.1, value=10.0, step=0.1)
    trim_s = st.number_input("先頭/末尾カット (秒)", min_value=0.0, value=5.0, step=0.5)
    unit_to_cm = st.number_input("単位→cm 係数 (mm→cm=0.1, cmなら1.0)", min_value=0.0001, value=0.1, step=0.1, format="%.4f")
    st.divider()
    st.subheader("図の表示範囲 (cm)")
    x_min, x_max = st.slider("ML (x) 範囲", -50.0, 50.0, (-6.0, 6.0), step=0.5)
    y_min, y_max = st.slider("AP (y) 範囲", -50.0, 50.0, (-6.0, 6.0), step=0.5)
    save_outputs = st.checkbox("outputs/ にファイル保存する", value=True)

cfg = CopConfig(fs=float(fs), lowpass_hz=float(lowpass), trim_seconds=float(trim_s), unit_to_cm=float(unit_to_cm))


st.markdown(
    """
**CSVの必要カラム**（大文字小文字は不問 / 前後のスペースはOK）  
- `time`（秒）  
- `COPx`（ML）  
- `COPy`（AP）
"""
)

uploaded = st.file_uploader("CSVをアップロード", type=["csv"])

# 既に app/data にCSVが置いてある場合は選択できるようにする
existing_csvs = sorted([p for p in DATA_DIR.glob("*.csv")])
selected_existing = None
if uploaded is None and existing_csvs:
    selected_existing = st.selectbox("または data/ 内のCSVから選択", ["(選択しない)"] + [p.name for p in existing_csvs])
    if selected_existing == "(選択しない)":
        selected_existing = None

csv_path: Path | None = None

if uploaded is not None:
    # アプリ内に保存
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    safe_name = uploaded.name.replace("/", "_").replace("\\", "_")
    csv_path = UPLOAD_DIR / f"{ts}_{safe_name}"
    csv_path.write_bytes(uploaded.getbuffer())
    st.success(f"アップロードCSVを保存しました: {csv_path.as_posix()}")
elif selected_existing:
    csv_path = DATA_DIR / selected_existing
    st.info(f"選択中: {csv_path.as_posix()}")

if csv_path is None:
    st.stop()

# ---- Analyze ----
try:
    res = analyze_cop_csv(csv_path, cfg=cfg)
except Exception as e:
    st.error(f"解析に失敗しました: {e}")
    st.stop()

metrics_df, bands_df = result_tables(res)

# ---- Outputs names ----
stem = csv_path.stem
ts = datetime.now().strftime("%Y%m%d_%H%M%S")
base = f"{stem}_{ts}"
panel_png_path = OUTPUT_DIR / f"{base}_cop_panel.png"
metrics_csv_path = OUTPUT_DIR / f"{base}_metrics.csv"
bands_csv_path = OUTPUT_DIR / f"{base}_bands.csv"

# ---- Layout ----
col1, col2 = st.columns([1, 1], gap="large")

with col1:
    st.subheader("数値（Metrics）")
    st.metric("95%信頼楕円 面積 (cm²)", f"{res.ellipse_area_cm2:.3f}")
    st.metric("外周面積（凸包） (cm²)", f"{res.hull_area_cm2:.3f}" if pd.notna(res.hull_area_cm2) else "nan")
    st.dataframe(metrics_df, use_container_width=True)

    st.subheader("帯域パワー比（0–3 Hz = 100%）")
    st.dataframe(bands_df, use_container_width=True)

with col2:
    st.subheader("可視化")
    fig = make_panel_figure(res, xlim=(x_min, x_max), ylim=(y_min, y_max))
    st.pyplot(fig, clear_figure=False)

    # figure bytes for download
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=300, bbox_inches="tight")
    buf.seek(0)

    st.download_button(
        label="図をPNGでダウンロード",
        data=buf,
        file_name=panel_png_path.name,
        mime="image/png",
    )

# ---- Save to app files ----
if save_outputs:
    # save png + csv to outputs/
    panel_png_path.write_bytes(buf.getvalue())
    metrics_df.to_csv(metrics_csv_path, index=False)
    bands_df.to_csv(bands_csv_path, index=False)

    st.success("outputs/ に保存しました")
    st.code(
        "\n".join(
            [
                panel_png_path.as_posix(),
                metrics_csv_path.as_posix(),
                bands_csv_path.as_posix(),
            ]
        )
    )

st.caption("※ Streamlit Community Cloud 等にデプロイした場合、ファイル保存は“実行環境内”で保持され、再デプロイ等で消えることがあります。永続保存が必要なら外部ストレージをご利用ください。")
