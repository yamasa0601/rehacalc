# BBS入力（脳卒中・大腿骨近位部骨折）

## できること
- BBS 14項目（各0〜4点、最大56点）をタップ入力
- 静的（1〜7）/動的（8〜14）/合計の自動計算
- 「シート1：脳卒中」「シート2：大腿骨近位部骨折」を選択
- CSVに入っている cut-off（閾値）を表示し、合計点で即判定
- 前回BBSを入れると差分Δを計算し、MDC/MCIDと照合
- 大腿骨近位部骨折シートはTUG・5回立ち座り（12秒以下）も任意で判定
- 入力は端末内に保存（PWA / オフライン可）

## GitHub Pagesで公開（おすすめ）
1. このフォルダ内の **index.html / manifest.json / service-worker.js / icon-*.png** をリポジトリ直下にアップロード
2. Settings → Pages → Deploy from a branch → main / (root)
3. 公開URL： https://<username>.github.io/<repo>/

## iPhoneで普段使い
Safariで公開URLを開く → 共有 → ホーム画面に追加
