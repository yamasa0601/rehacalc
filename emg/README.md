# EMG歩行解析（スマホ対応・サーバ不要）

`rehacalc` の GitHub Pages / PWA にそのまま追加できる、ブラウザ完結のEMG解析ページです。

## できること
- CSVを2つ読み込み（左右など）
- 50–450Hz bandpass（HP→LP、2次Butterworth×2、forward-backwardで位相補正）
- 整流 → 50ms RMS包絡
- バースト立ち上がりからHS推定（median + k×MAD、kは歩行らしさで自動選択）
- %peak（RMS最大値=100%）
- HS–HSで周期分割 → 0–100%に正規化（N点）
- 平均±SDプロット、PNG保存、HS一覧CSV保存

## 配置
リポジトリ直下に `emg/` フォルダとして置いてください。

例：
```
rehacalc/
  index.html
  service-worker.js
  manifest.json
  emg/
    index.html
    app.js
    style.css
```

## GitHub Pages で開く
公開URLが `https://<user>.github.io/rehacalc/` の場合：
- `https://<user>.github.io/rehacalc/emg/`

## PWAとしてオフライン運用したい場合
既存の `service-worker.js` がプリキャッシュ方式なら、キャッシュ対象に以下を追加してください：
- `emg/index.html`
- `emg/app.js`
- `emg/style.css`

PlotlyはCDN読み込みです。完全オフライン化したい場合はPlotlyをリポジトリに同梱してください。
