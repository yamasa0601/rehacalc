
# EMG Heel-Strike Estimator (Walking Only)

2つのヘッダー付きEMG CSVをアップロードすると、以下を自動で実行してプロットまで出します。

- 50–450 Hz bandpass
- 整流
- 50 ms RMS（包絡）
- **%peak 正規化**（RMS最大値 = 100%）
- **筋活動バーストの立ち上がり**から Heel Strike（HS）推定（median + k×MAD の動的閾値）
- HS間で歩行周期を切り出し、0–100%に時間正規化して平均±SDを表示
- HS一覧CSVのダウンロード / 図（PNG）の保存

> 想定：歩行のみ。サンプリング周波数はデフォルト1000 Hz（サイドバーで変更可能）。

---

## Install

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
# macOS/Linux:
source .venv/bin/activate

pip install -r requirements.txt
```

## Run

```bash
streamlit run app.py
```

ブラウザが開くので、CSVを2つアップロードして「解析してプロット」を押してください。  
結果は画面表示され、必要なら `outputs/` にPNGとして保存されます。

---

## CSV format

このアプリは2種類を自動判別します：

### A) センサ出力CSV（先頭にメタ情報、途中から `No,電圧 mV,event` が始まる）
あなたの例（sensor1_analog.csv 等）に対応しています。

### B) 普通の表形式CSV
例：`time, value` のような列。時間列が無ければ先頭列をサンプル/インデックスとして扱います。

---

## HS推定はTA以外でもできる？

**できることが多い**ですが、**必ず成功するわけではありません**。

- TA：遊脚終期〜接地付近でバーストが出やすく、HS推定に向くことが多い
- MG/腓腹筋：初期立脚〜推進期で活動が強く、HSより後ろ側のイベントに寄りがち
- RF/BF：速度や歩容でパターンが変わりやすい

本アプリは「バーストの立ち上がり」をHSとして扱うため、筋の種類や貼付位置でズレが出ます。  
ズレる場合は **`HSオフセット(ms)`** を調整してください。  
（本当にHSが必要なら、フットスイッチ/IMU/床反力などの併用が最も確実です。）

---

## License

MIT
