# COP analysis (Java core)

Pythonノートブックの処理をJavaで再現するための最小構成です。  
CLIとして動かせるので、ここを土台にAndroid/JavaFX/サーバーアプリへ展開できます。

## 使い方

```
mvn -q package
java -jar target/cop-java-1.0.0.jar /path/to/OE.csv
```

## 既定パラメータ

- FS=100 Hz
- LPF=10 Hz (4次Butterworth, filtfilt)
- 前後除外=5 秒
- 入力単位=mm, 出力単位=cm
- Welch: Hann, 50% overlap, detrend=mean

必要なら `Main.java` のパラメータを変更してください。

## GitHub用メモ

- `cop-java/` をそのままリポジトリに置けばOK
- GitHub Actionsでビルドする場合は `mvn -q package`



