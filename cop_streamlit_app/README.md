# COP Streamlit App

CSVをアップロード（または `data/` に置くだけ）で、COPの可視化と数値を自動計算するStreamlitアプリです。

## できること
- 先頭/末尾の指定秒数をトリミング
- COPを平均0にセンタリング
- ローパスフィルタ（デフォルト10Hz）
- 単位をcmへ変換（mm→cmなら0.1、すでにcmなら1.0）
- 95%信頼楕円面積（cm²）
- 凸包（外周）面積（cm²）
- 帯域パワー比（0–3Hzの総パワー=100%）
- 図と数値（CSV/PNG）を `outputs/` に保存

## CSVフォーマット
必要カラム（大文字小文字は不問、前後スペースはOK）:
- `time` : 秒
- `COPx` : ML
- `COPy` : AP

## ローカル実行
```bash
python -m venv .venv
# mac/linux
source .venv/bin/activate
# windows: .venv\Scripts\activate

pip install -r requirements.txt
streamlit run app.py
```

## データを置く場所
- アップロードしたCSVは `data/uploads/` に保存されます
- 手動で置くなら `data/` にCSVを置いてください（アプリ上で選択できます）
- 解析結果は `outputs/` に保存されます

## GitHubで使う（例）
1. このフォルダ一式をGitHubにpush
2. ローカルなら上記コマンドで実行  
3. Streamlit Community Cloud等にデプロイする場合は、リポジトリ内の `app.py` をエントリに指定してください
