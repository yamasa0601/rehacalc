import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

public class CopAnalyzer {
    public static final class Params {
        public double fs = 100.0;
        public double cutoffHz = 10.0;
        public double trimSec = 5.0;
        public String unitIn = "mm";   // "mm" or "cm"
        public String unitOut = "cm";  // "mm" or "cm"
    }

    public static final class Result {
        public double ellipseArea;
        public double hullArea;
        public double lfRatio;
        public double mfRatio;
        public double hfRatio;
        public int nPoints;
        public double durationSec;
        public double[] x;
        public double[] y;
    }

    public static Result analyzeCsv(Path csvPath, Params params) throws IOException {
        CsvData data = readCsv(csvPath);
        return analyze(data.time, data.copx, data.copy, params);
    }

    public static Result analyze(double[] time, double[] copx, double[] copy, Params params) {
        if (time.length < 2) {
            throw new IllegalArgumentException("データ点が少なすぎます。");
        }

        // 前後除外（秒）
        double t0 = time[0] + params.trimSec;
        double t1 = time[time.length - 1] - params.trimSec;
        Trimmed trimmed = trimByTime(time, copx, copy, t0, t1);

        if (trimmed.time.length < 20) {
            throw new IllegalArgumentException("前後除外後にデータが少なすぎます。");
        }

        // センタリング（生COP）
        double[] cx = center(trimmed.copx);
        double[] cy = center(trimmed.copy);

        // 4次Butterworth ローパス（filtfilt）
        double[] fx = filtfiltButter4Lowpass(cx, params.fs, params.cutoffHz);
        double[] fy = filtfiltButter4Lowpass(cy, params.fs, params.cutoffHz);

        // 単位変換（Pythonと同じ：ローパス後に変換）
        double unitToCm = "mm".equalsIgnoreCase(params.unitIn) ? 0.1 : 1.0;
        double cmToOut = "mm".equalsIgnoreCase(params.unitOut) ? 10.0 : 1.0;
        double unitToOut = unitToCm * cmToOut;
        double[] xOut = scale(fx, unitToOut);
        double[] yOut = scale(fy, unitToOut);

        // 95%楕円面積
        double[][] cov = cov2(xOut, yOut);
        Ellipse ell = ellipseFromCov(cov);
        double ellipseArea = ell.area;

        // 凸包面積
        double[][] hull = convexHull(xOut, yOut);
        double hullArea = polygonArea(hull);

        // Welch PSD（0–3Hz）
        WelchPSD psdX = welchPSD(fx, params.fs);
        WelchPSD psdY = welchPSD(fy, params.fs);
        double[] f = psdX.f;
        double[] pTotal = new double[f.length];
        for (int i = 0; i < f.length; i++) pTotal[i] = psdX.p[i] + psdY.p[i];

        BandRatios ratios = bandRatios(f, pTotal);

        Result res = new Result();
        res.ellipseArea = ellipseArea;
        res.hullArea = hullArea;
        res.lfRatio = ratios.lf;
        res.mfRatio = ratios.mf;
        res.hfRatio = ratios.hf;
        res.nPoints = trimmed.time.length;
        res.durationSec = trimmed.time[trimmed.time.length - 1] - trimmed.time[0];
        res.x = xOut;
        res.y = yOut;
        return res;
    }

    // ===== CSV =====
    private static final class CsvData {
        double[] time;
        double[] copx;
        double[] copy;
    }

    private static CsvData readCsv(Path path) throws IOException {
        List<String> lines = Files.readAllLines(path, StandardCharsets.UTF_8);
        if (lines.isEmpty()) throw new IllegalArgumentException("CSVが空です。");

        String[] headers = splitCsvLine(lines.get(0));
        int timeIdx = findHeader(headers, "time", "t", "timestamp", "sec", "seconds");
        int xIdx = findHeader(headers, "copx", "cop_x", "cop-x", "ml", "x");
        int yIdx = findHeader(headers, "copy", "copyy", "cop_y", "cop-y", "ap", "y");

        if (timeIdx < 0 || xIdx < 0 || yIdx < 0) {
            throw new IllegalArgumentException("必須列が見つかりません。");
        }

        List<Double> time = new ArrayList<>();
        List<Double> copx = new ArrayList<>();
        List<Double> copy = new ArrayList<>();

        for (int i = 1; i < lines.size(); i++) {
            String line = lines.get(i).trim();
            if (line.isEmpty()) continue;
            String[] cols = splitCsvLine(line);
            if (cols.length <= Math.max(timeIdx, Math.max(xIdx, yIdx))) continue;

            double t = toNumber(cols[timeIdx]);
            double x = toNumber(cols[xIdx]);
            double y = toNumber(cols[yIdx]);
            if (!Double.isFinite(t) || !Double.isFinite(x) || !Double.isFinite(y)) continue;

            time.add(t);
            copx.add(x);
            copy.add(y);
        }

        if (time.size() < 2) throw new IllegalArgumentException("数値として読める行がありません。");

        CsvData data = new CsvData();
        data.time = toArray(time);
        data.copx = toArray(copx);
        data.copy = toArray(copy);
        return data;
    }

    private static String[] splitCsvLine(String line) {
        // ここでは簡易分割（ダブルクォートなしのCSV想定）
        String[] raw = line.split(",");
        for (int i = 0; i < raw.length; i++) raw[i] = raw[i].trim();
        return raw;
    }

    private static int findHeader(String[] headers, String... candidates) {
        for (int i = 0; i < headers.length; i++) {
            String h = normalizeKey(headers[i]);
            for (String c : candidates) {
                if (h.equals(normalizeKey(c))) return i;
            }
        }
        return -1;
    }

    private static String normalizeKey(String s) {
        return s == null ? "" : s.trim().toLowerCase().replaceAll("\\s+", "");
    }

    private static double toNumber(String v) {
        if (v == null) return Double.NaN;
        String s0 = v.trim();
        if (s0.isEmpty()) return Double.NaN;
        String s = s0;
        if (s.contains(",") && !s.contains(".")) {
            String[] parts = s.split(",");
            if (parts.length == 2 && parts[1].length() <= 3) s = parts[0] + "." + parts[1];
        }
        String m = s.replaceAll("^.*?([-+]?\\d*\\.?\\d+(?:[eE][-+]?\\d+)?).*$", "$1");
        try {
            double num = Double.parseDouble(m);
            return Double.isFinite(num) ? num : Double.NaN;
        } catch (NumberFormatException e) {
            return Double.NaN;
        }
    }

    private static double[] toArray(List<Double> list) {
        double[] arr = new double[list.size()];
        for (int i = 0; i < list.size(); i++) arr[i] = list.get(i);
        return arr;
    }

    // ===== preprocessing =====
    private static final class Trimmed {
        double[] time;
        double[] copx;
        double[] copy;
    }

    private static Trimmed trimByTime(double[] time, double[] x, double[] y, double t0, double t1) {
        List<Double> nt = new ArrayList<>();
        List<Double> nx = new ArrayList<>();
        List<Double> ny = new ArrayList<>();
        for (int i = 0; i < time.length; i++) {
            double t = time[i];
            if (t >= t0 && t <= t1) {
                nt.add(t);
                nx.add(x[i]);
                ny.add(y[i]);
            }
        }
        Trimmed tr = new Trimmed();
        tr.time = toArray(nt);
        tr.copx = toArray(nx);
        tr.copy = toArray(ny);
        return tr;
    }

    private static double[] center(double[] arr) {
        double m = mean(arr);
        double[] out = new double[arr.length];
        for (int i = 0; i < arr.length; i++) out[i] = arr[i] - m;
        return out;
    }

    private static double mean(double[] arr) {
        double s = 0;
        for (double v : arr) s += v;
        return s / arr.length;
    }

    private static double[] scale(double[] arr, double k) {
        double[] out = new double[arr.length];
        for (int i = 0; i < arr.length; i++) out[i] = arr[i] * k;
        return out;
    }

    // ===== Butterworth 4th-order filtfilt =====
    private static final double[] B_FIXED = {0.004824343357716228, 0.01929737343086491, 0.028946060146297367, 0.01929737343086491, 0.004824343357716228};
    private static final double[] A_FIXED = {1.0, -2.369513007182037, 2.3139884144158806, -1.0546654058785673, 0.18737949236818535};
    private static final double[] ZI_FIXED = {0.9951756566422838, -1.3936347201975356, 0.891407628579869, -0.1825551467889937};

    private static double[] filtfiltButter4Lowpass(double[] x, double fs, double cutoff) {
        if (x.length <= 15) return Arrays.copyOf(x, x.length);
        double wn = cutoff / (0.5 * fs);
        if (Math.abs(wn - 0.2) <= 1e-3) {
            return filtfiltIir(x, B_FIXED, A_FIXED, ZI_FIXED);
        }
        Biquad[] biquads = designButterworth4LowpassBiquads(fs, cutoff);
        double[] y = Arrays.copyOf(x, x.length);
        for (Biquad bq : biquads) y = filtfiltBiquad(y, bq);
        return y;
    }

    private static double[] filtfiltIir(double[] x, double[] b, double[] a, double[] ziBase) {
        int padlen = 3 * (Math.max(a.length, b.length) - 1);
        if (x.length <= padlen) return Arrays.copyOf(x, x.length);
        double[] ext = oddExtend(x, padlen);
        double[] y1 = lfilter(ext, b, a, ziBase, ext[0]);
        double[] rev = reverse(y1);
        double[] y2 = lfilter(rev, b, a, ziBase, rev[0]);
        double[] y2r = reverse(y2);
        return Arrays.copyOfRange(y2r, padlen, padlen + x.length);
    }

    private static double[] oddExtend(double[] x, int padlen) {
        int n = x.length;
        double[] left = new double[padlen];
        double[] right = new double[padlen];
        for (int i = 1; i <= padlen; i++) left[i - 1] = 2 * x[0] - x[i];
        for (int i = 1; i <= padlen; i++) right[i - 1] = 2 * x[n - 1] - x[n - 1 - i];
        reverseInPlace(left);
        reverseInPlace(right);
        double[] out = new double[padlen + n + padlen];
        System.arraycopy(left, 0, out, 0, padlen);
        System.arraycopy(x, 0, out, padlen, n);
        System.arraycopy(right, 0, out, padlen + n, padlen);
        return out;
    }

    private static double[] lfilter(double[] x, double[] b, double[] a, double[] ziBase, double x0) {
        int order = Math.max(b.length, a.length) - 1;
        double[] z = new double[order];
        for (int i = 0; i < order; i++) z[i] = ziBase[i] * x0;

        double[] y = new double[x.length];
        for (int n = 0; n < x.length; n++) {
            double xn = x[n];
            double yn = b[0] * xn + z[0];
            for (int i = 1; i < order; i++) {
                double bi = i < b.length ? b[i] : 0;
                double ai = i < a.length ? a[i] : 0;
                z[i - 1] = bi * xn + z[i] - ai * yn;
            }
            double bLast = order < b.length ? b[order] : 0;
            double aLast = order < a.length ? a[order] : 0;
            z[order - 1] = bLast * xn - aLast * yn;
            y[n] = yn;
        }
        return y;
    }

    private static double[] reverse(double[] arr) {
        double[] out = Arrays.copyOf(arr, arr.length);
        reverseInPlace(out);
        return out;
    }

    private static void reverseInPlace(double[] arr) {
        for (int i = 0, j = arr.length - 1; i < j; i++, j--) {
            double t = arr[i];
            arr[i] = arr[j];
            arr[j] = t;
        }
    }

    private static final class Biquad {
        double b0, b1, b2, a1, a2;
    }

    private static Biquad[] designButterworth4LowpassBiquads(double fs, double fc) {
        double[] Qs = {0.5411961, 1.3065630};
        Biquad[] bqs = new Biquad[Qs.length];
        for (int i = 0; i < Qs.length; i++) bqs[i] = designRBJLowpass(fs, fc, Qs[i]);
        return bqs;
    }

    private static Biquad designRBJLowpass(double fs, double fc, double Q) {
        double w0 = 2 * Math.PI * fc / fs;
        double cosw0 = Math.cos(w0);
        double sinw0 = Math.sin(w0);
        double alpha = sinw0 / (2 * Q);
        double b0 = (1 - cosw0) / 2;
        double b1 = 1 - cosw0;
        double b2 = (1 - cosw0) / 2;
        double a0 = 1 + alpha;
        double a1 = -2 * cosw0;
        double a2 = 1 - alpha;

        b0 /= a0; b1 /= a0; b2 /= a0; a1 /= a0; a2 /= a0;
        Biquad bq = new Biquad();
        bq.b0 = b0; bq.b1 = b1; bq.b2 = b2; bq.a1 = a1; bq.a2 = a2;
        return bq;
    }

    private static double[] filtfiltBiquad(double[] x, Biquad bq) {
        int padlen = 12;
        if (x.length <= padlen) return Arrays.copyOf(x, x.length);
        double[] ext = oddExtend(x, padlen);
        double[] y1 = biquadFilter(ext, bq);
        double[] y2 = biquadFilter(reverse(y1), bq);
        double[] y2r = reverse(y2);
        return Arrays.copyOfRange(y2r, padlen, padlen + x.length);
    }

    private static double[] biquadFilter(double[] x, Biquad bq) {
        double[] y = new double[x.length];
        double z1 = 0, z2 = 0;
        for (int i = 0; i < x.length; i++) {
            double xn = x[i];
            double yn = bq.b0 * xn + z1;
            z1 = bq.b1 * xn + z2 - bq.a1 * yn;
            z2 = bq.b2 * xn - bq.a2 * yn;
            y[i] = yn;
        }
        return y;
    }

    // ===== covariance / ellipse =====
    private static double[][] cov2(double[] x, double[] y) {
        int n = x.length;
        double mx = mean(x), my = mean(y);
        double sxx = 0, syy = 0, sxy = 0;
        for (int i = 0; i < n; i++) {
            double dx = x[i] - mx;
            double dy = y[i] - my;
            sxx += dx * dx;
            syy += dy * dy;
            sxy += dx * dy;
        }
        double denom = n - 1;
        return new double[][]{{sxx / denom, sxy / denom}, {sxy / denom, syy / denom}};
    }

    private static final class Ellipse {
        double area;
    }

    private static Ellipse ellipseFromCov(double[][] cov) {
        double a = cov[0][0], b = cov[0][1], d = cov[1][1];
        double tr = a + d;
        double disc = Math.sqrt((a - d) * (a - d) + 4 * b * b);
        double l1 = (tr + disc) / 2;
        double l2 = (tr - disc) / 2;
        double chi2ppf = 5.991464547107979;
        double area = Math.PI * chi2ppf * Math.sqrt(Math.max(l1, 0) * Math.max(l2, 0));
        Ellipse e = new Ellipse();
        e.area = area;
        return e;
    }

    // ===== convex hull =====
    private static double[][] convexHull(double[] x, double[] y) {
        int n = x.length;
        double[][] pts = new double[n][2];
        for (int i = 0; i < n; i++) {
            pts[i][0] = x[i];
            pts[i][1] = y[i];
        }
        Arrays.sort(pts, Comparator.<double[]>comparingDouble(p -> p[0]).thenComparingDouble(p -> p[1]));
        List<double[]> lower = new ArrayList<>();
        for (double[] p : pts) {
            while (lower.size() >= 2 && cross(lower.get(lower.size() - 2), lower.get(lower.size() - 1), p) <= 0) {
                lower.remove(lower.size() - 1);
            }
            lower.add(p);
        }
        List<double[]> upper = new ArrayList<>();
        for (int i = pts.length - 1; i >= 0; i--) {
            double[] p = pts[i];
            while (upper.size() >= 2 && cross(upper.get(upper.size() - 2), upper.get(upper.size() - 1), p) <= 0) {
                upper.remove(upper.size() - 1);
            }
            upper.add(p);
        }
        if (!lower.isEmpty()) lower.remove(lower.size() - 1);
        if (!upper.isEmpty()) upper.remove(upper.size() - 1);
        List<double[]> hull = new ArrayList<>(lower);
        hull.addAll(upper);
        return hull.toArray(new double[0][0]);
    }

    private static double cross(double[] o, double[] a, double[] b) {
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]);
    }

    private static double polygonArea(double[][] poly) {
        if (poly == null || poly.length < 3) return 0;
        double s = 0;
        for (int i = 0; i < poly.length; i++) {
            double[] p1 = poly[i];
            double[] p2 = poly[(i + 1) % poly.length];
            s += p1[0] * p2[1] - p2[0] * p1[1];
        }
        return Math.abs(s) / 2;
    }

    // ===== Welch PSD =====
    private static final class WelchPSD {
        double[] f;
        double[] p;
    }

    private static WelchPSD welchPSD(double[] x, double fs) {
        int N = x.length;
        int nperseg = Math.min(1024, N);
        int noverlap = nperseg / 2;
        int step = nperseg - noverlap;
        if (step <= 0) throw new IllegalArgumentException("noverlap too large");

        double[] w = hann(nperseg);
        double w2sum = 0;
        for (double v : w) w2sum += v * v;

        FFT fft = new FFT(nperseg);
        double[] out = fft.createComplexArray();
        double[] input = fft.createComplexArray();

        int nfreq = nperseg / 2 + 1;
        double[] P = new double[nfreq];
        int nseg = 0;

        for (int start = 0; start + nperseg <= N; start += step) {
            double m = 0;
            for (int i = 0; i < nperseg; i++) m += x[start + i];
            m /= nperseg;

            for (int i = 0; i < nperseg; i++) {
                double v = (x[start + i] - m) * w[i];
                input[2 * i] = v;
                input[2 * i + 1] = 0;
            }

            fft.transform(out, input);
            double scale = 1 / (fs * w2sum);
            for (int k = 0; k < nfreq; k++) {
                double re = out[2 * k];
                double im = out[2 * k + 1];
                double pk = (re * re + im * im) * scale;
                if (k != 0 && k != nperseg / 2) pk *= 2;
                P[k] += pk;
            }
            nseg += 1;
        }
        if (nseg == 0) throw new IllegalArgumentException("データが短すぎてWelch計算できません。");
        for (int k = 0; k < P.length; k++) P[k] /= nseg;

        double[] f = new double[nfreq];
        for (int k = 0; k < nfreq; k++) f[k] = (k * fs) / nperseg;

        WelchPSD psd = new WelchPSD();
        psd.f = f;
        psd.p = P;
        return psd;
    }

    private static double[] hann(int N) {
        double[] w = new double[N];
        for (int n = 0; n < N; n++) w[n] = 0.5 - 0.5 * Math.cos(2 * Math.PI * n / (N - 1));
        return w;
    }

    private static final class BandRatios {
        double lf;
        double mf;
        double hf;
    }

    private static BandRatios bandRatios(double[] f, double[] p) {
        BandRatios r = new BandRatios();
        r.lf = bandRatio(f, p, 0.0, 0.3);
        r.mf = bandRatio(f, p, 0.3, 1.0);
        r.hf = bandRatio(f, p, 1.0, 3.0);
        return r;
    }

    private static double bandRatio(double[] f, double[] p, double a, double b) {
        double total = trapz(f, p, 0.0, 3.0);
        if (!(total > 0)) return Double.NaN;
        double band = trapz(f, p, a, b);
        return band / total * 100.0;
    }

    private static double trapz(double[] f, double[] p, double fmin, double fmax) {
        double s = 0;
        for (int i = 0; i < f.length - 1; i++) {
            if (f[i] < fmin || f[i + 1] > fmax) continue;
            double dx = f[i + 1] - f[i];
            s += dx * (p[i] + p[i + 1]) / 2;
        }
        return s;
    }

    // ===== FFT (Cooley-Tukey, pow2 only; DFT fallback) =====
    private static final class FFT {
        private final int size;
        private final boolean isPow2;
        private final double[] cos;
        private final double[] sin;

        FFT(int size) {
            this.size = size;
            this.isPow2 = (size & (size - 1)) == 0 && size > 0;
            if (isPow2) {
                cos = new double[size / 2];
                sin = new double[size / 2];
                for (int i = 0; i < size / 2; i++) {
                    double ang = -2 * Math.PI * i / size;
                    cos[i] = Math.cos(ang);
                    sin[i] = Math.sin(ang);
                }
            } else {
                cos = null;
                sin = null;
            }
        }

        double[] createComplexArray() {
            return new double[size * 2];
        }

        void transform(double[] out, double[] input) {
            if (isPow2) {
                fft(out, input);
            } else {
                dft(out, input);
            }
        }

        private void fft(double[] out, double[] input) {
            int n = size;
            for (int i = 0; i < n * 2; i++) out[i] = input[i];

            int j = 0;
            for (int i = 0; i < n; i++) {
                if (i < j) {
                    int ia = 2 * i, ja = 2 * j;
                    double tr = out[ia], ti = out[ia + 1];
                    out[ia] = out[ja]; out[ia + 1] = out[ja + 1];
                    out[ja] = tr; out[ja + 1] = ti;
                }
                int m = n >> 1;
                while (j >= m && m >= 2) { j -= m; m >>= 1; }
                j += m;
            }

            for (int size = 2; size <= n; size <<= 1) {
                int half = size >> 1;
                int step = n / size;
                for (int i = 0; i < n; i += size) {
                    for (int k = 0; k < half; k++) {
                        int idx = k * step;
                        double c = cos[idx];
                        double s = sin[idx];
                        int a = 2 * (i + k);
                        int b = 2 * (i + k + half);
                        double br = out[b];
                        double bi = out[b + 1];
                        double tr = c * br - s * bi;
                        double ti = s * br + c * bi;
                        double ur = out[a];
                        double ui = out[a + 1];
                        out[a] = ur + tr;
                        out[a + 1] = ui + ti;
                        out[b] = ur - tr;
                        out[b + 1] = ui - ti;
                    }
                }
            }
        }

        private void dft(double[] out, double[] input) {
            int n = size;
            for (int k = 0; k < n; k++) {
                double re = 0;
                double im = 0;
                for (int t = 0; t < n; t++) {
                    double xr = input[2 * t];
                    double xi = input[2 * t + 1];
                    double ang = -2 * Math.PI * k * t / n;
                    double c = Math.cos(ang);
                    double s = Math.sin(ang);
                    re += xr * c - xi * s;
                    im += xr * s + xi * c;
                }
                out[2 * k] = re;
                out[2 * k + 1] = im;
            }
        }
    }
}



