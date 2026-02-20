import java.nio.file.Path;

public class Main {
    public static void main(String[] args) throws Exception {
        if (args.length < 1) {
            System.out.println("Usage: java Main /path/to/OE.csv");
            return;
        }

        CopAnalyzer.Params params = new CopAnalyzer.Params();
        params.fs = 100.0;
        params.cutoffHz = 10.0;
        params.trimSec = 5.0;
        params.unitIn = "mm";   // CSVがmmならmm、cmならcm
        params.unitOut = "cm";  // 出力をcmにする

        CopAnalyzer.Result res = CopAnalyzer.analyzeCsv(Path.of(args[0]), params);

        System.out.printf("95%% confidence ellipse area: %.3f %s^2%n",
                res.ellipseArea, params.unitOut);
        System.out.printf("Outer area (convex hull): %.3f %s^2%n",
                res.hullArea, params.unitOut);
        System.out.printf("Band power ratios (0–3 Hz total = 100%%):%n");
        System.out.printf("  LF_0-0.3: %.2f %%\n", res.lfRatio);
        System.out.printf("  MF_0.3-1.0: %.2f %%\n", res.mfRatio);
        System.out.printf("  HF_1-3: %.2f %%\n", res.hfRatio);
    }
}



