import java.util.Objects;
import java.util.Random;
import java.util.Scanner;
import java.util.concurrent.TimeUnit;
import java.util.function.Function;

public class Main {

    public static void main(String[] args) throws InterruptedException {

        //Float();

        //Double();

        // xxx();


        //zzz();

        /*bisekcja(-4,-3);

        System.out.println();

        bisekcja(-2,-1);

        System.out.println();

        bisekcja(-0,2);

        System.out.println();

        bisekcja(4,5);*/

        //zadanie1();

        //zadanie2();

        /*przykład1();

        przykład2();

        przykład3();

        przykład4();*/

        //test();

        //monte();

        //monteBlad();

        //pi();

        //System.out.println(Math.atan(1)*4);

        //x();

        //randomWalker2D(10);

        //y();

        //licz();

        //OHzad2(0, 20 * Math.PI, 10);

        //OH(0,1,10);

        //Verlet();

        lab8();
    }

    public static void test() {

        double a = 0.0;

        double b = 1.0;

        int n = 10;

        for (int i = 0; i < 9; i++) {

            double result1 = trapezoidalMethod(a, b, n);

            System.out.println(result1);
            n = n * 10;
        }

/// policzenia wyniku całkowania dla każdej z funkcji testowych

        double result = trapezoidalMethod(a, b, n);

    }

    public static double trapezoidalMethod(double a, double b, int n) {// Eugodt poprzeszkolu

        double h = (b - a) / n;

        double sum = (zad2(a) + zad2(b)) / 2.0; // suna wartości funkcji na krancach przedziału

        for (int i = 1; i < n; i++) {
            double x = a + i * h; // punkt środkowy 1-tegu pooprzedziału

            sum += zad2(x); // dodanie wartotel funkcji w punkcie środkowye de suny

        }

        return h * sum;
        // wynik satuwaele przybliżony, autody Trapezow


    }


    public static void xxx() {
        double[] x = {2.436,
                2.435,
                2.442,
                2.448,
                2.445,
                2.439,
                2.445,
                2.431,
                2.435,
                2.442
        };

        double wynik = 0;

        for (int i = 0; i < x.length; i++) {
            wynik += Math.pow(x[i] - 2.44, 2);
        }

        System.out.format("%f", wynik);


    }


    public static void zzz() {
        double[] x = {2.433,
                2.475,
                2.435,
                2.451,
                2.441,
                2.445,
                2.439,
                2.44,
                2.445,
                2.518,
        };

        double wynik = 0;

        for (int i = 0; i < x.length; i++) {
            wynik += Math.pow(x[i] - 2.45, 2);
        }

        System.out.format("%f", wynik);


    }

    public static void Float() {
        long millisActualTime = System.currentTimeMillis();

        float a = 2;

        float b = 3;

        float c = 6;

        float d = 7;

        float x;

        float y;

        float z;

        float t;

        float u;

        double time = 1 * Math.pow(10, 9);

        for (int i = 0; i < time; i++) {
            x = a / c;

            y = a / b;

            z = d / c;

            t = (b / d) * c;

            u = (b * c) / d;

            a++;

            b++;

            c++;

            d++;

        }


        long executionTime = System.currentTimeMillis() - millisActualTime;

        System.out.println((double) executionTime / 1000);
    }

    public static void Double() {
        long millisActualTime = System.currentTimeMillis();

        double a = 2;

        double b = 3;

        double c = 6;

        double d = 7;

        double x;

        double y;

        double z;

        double t;

        double u;

        double time = 1 * Math.pow(10, 9);

        for (int i = 0; i < time; i++) {
            x = a / c;

            y = a / b;

            z = d / c;

            t = (b / d) * c;

            u = (b * c) / d;

            a++;

            b++;

            c++;

            d++;

        }


        long executionTime = System.currentTimeMillis() - millisActualTime;

        System.out.println((double) executionTime / 1000);
    }

    public static void bisekcja(double a, double b) {

        //x^4-15x^2-15x+5;

        double eps = 0.00000000000001;

        double c = 0;

        double funA, funB, funC;

        funA = funkcja(a);

        funB = funkcja(b);

        if (funA * funB >= 0) {
            throw new RuntimeException("funkcja nie ma miejsc zerowych");
        }

        while (Math.abs(b - a) > eps) {

            c = (a + b) / 2;

            funC = funkcja(c);

            if (funA * funC < 0) {
                b = c;
            } else {
                a = c;
                funA = funC;
            }

        }

        System.out.print("Miejsce zerowe to ");

        System.out.format("%.16f", c);

    }

    public static double funkcja(double wartość) {

        double v = Math.pow(wartość, 4) - 15 * Math.pow(wartość, 2) - 15 * wartość + 5;
        return v;
    }

    public static void przykład1() {
        System.out.println(metodaTrapezów(0, 10, 0));

        System.out.println(metodaTrapezów(-1, 1, 1));

        System.out.println(metodaTrapezów(-10, 10, 2));

        System.out.println(metodaTrapezów(-10, 0, 3));

        System.out.println(metodaTrapezów(0, 1, 4));

    }

    public static void przykład2() {
        System.out.println(metodaProstokątów(0, 10, 0));

        System.out.println(metodaProstokątów(-1, 1, 1));

        System.out.println(metodaProstokątów(-10, 10, 2));

        System.out.println(metodaProstokątów(-10, 0, 3));

        System.out.println(metodaProstokątów(0, 1, 4));

    }

    public static void przykład3() {
        System.out.println(Simpson(0, 10, 0));

        System.out.println(Simpson(-1, 1, 1));

        System.out.println(Simpson(-10, 10, 2));

        System.out.println(Simpson(-10, 0, 3));

        System.out.println(Simpson(0, 1, 4));

    }

    public static void przykład4() {

        System.out.println("wynik całek    wynik metody Gaus");

        double aT = Gaus(0, 10, 0);

        System.out.println("90              " + aT);

        double bT = Gaus(-1, 1, 1);

        System.out.println("-2              " + bT);

        double cT = Gaus(-10, 10, 2);

        System.out.println("686.666         " + cT);

        double dT = Gaus(-10, 0, 3);

        System.out.println("-10 100         " + dT);

        double eT = Gaus(0, 1, 4);

        System.out.printf("0               %,.32f", eT);
        System.out.println();

        double fT = Gaus(0, 1, 5);

        System.out.printf("0               %,.32f", fT);
        System.out.println();

        double gT = Gaus(0, 1, 6);

        System.out.printf("0               %,.17f", gT);
        System.out.println();

    }

    public static double metodaTrapezów(double xp, double xk, int pow) {
        double calka = xk - xp;

        calka *= (func(xp, pow) + func(xk, pow));

        calka /= 2;

        return calka;
    }

    public static double metodaProstokątów(double xp, double xk, int pow) {

        double calka = xk - xp;

        double x = xp + xk;

        calka *= func(x / 2, pow);

        return calka;
    }

    public static double Simpson(double xp, double xk, int pow) {

        double calka = (xk - xp) / 6;

        double x = xp + xk;

        double calka2 = func(xp, pow) + 4 * func(x / 2, pow) + func(xk, pow);

        calka *= calka2;

        return calka;
    }

    public static double Gaus(double xp, double xk, int pow) {
        double[] x = {-0.906179845938663992797626878299392965125651910762530862873762286,
                -0.538469310105683091036314420700208804967286606905559956202231627,
                0,
                0.538469310105683091036314420700208804967286606905559956202231627,
                0.906179845938663992797626878299392965125651910762530862873762286};
        double[] w = {0.2369268850561890875142640407199173626432600022124140155828278882
                , 0.4786286704993664680412915148356381929122955533431415399727276673
                , 0.5688888888888888888888888888888888888888888888888888888888888888
                , 0.4786286704993664680412915148356381929122955533431415399727276673
                , 0.2369268850561890875142640407199173626432600022124140155828278882};

        double calka = 0;

        double temp;

        double temp2;

        for (int i = 0; i < 5; i++) {

            temp = (xk - xp) / 2;

            temp *= x[i];

            temp2 = (xk + xp) / 2;

            temp = func(temp + temp2, pow);

            temp *= w[i];

            calka += temp;
        }

        calka *= ((xk - xp) / 2);

        return calka;
    }

    public static double func(double x, int pow) {

        if (pow == 0) {
            return fA(x);
        } else if (pow == 1) {
            return fB(x);
        } else if (pow == 2) {
            return fC(x);
        } else if (pow == 3) {
            return fD(x);
        } else if (pow == 4) {
            return fE(x);
        } else if (pow == 5) {
            return fF(x);
        } else if (pow == 6) {
            return fG(x);
        }

        return 0;
    }

    public static double fA(double x) {
        return 9;
    }

    public static double fB(double x) {
        return (3 * x) - 1;
    }

    public static double fC(double x) {
        return Math.pow(x, 2) + 1;
    }

    public static double fD(double x) {
        return 4 * Math.pow(x, 3) + 2 * x;
    }

    public static double fE(double x) {
        return 5 * Math.pow(x, 4) - 2 * x;
    }

    public static double fF(double x) {
        return 10 * Math.pow(x, 9) - 8 * Math.pow(x, 7);
    }

    public static double fG(double x) {
        return 15 * Math.pow(x, 14) - 12 * Math.pow(x, 11);
    }

    public static double zad2(double x) {
        return 354 * Math.pow(x, 117);
    }

    private static void zadanie2() {

        System.out.println("  N               MT                     MP                       MS                     MGL");

        double t10 = MT(10, 1, 0);

        double p10 = MP(10, 1, 0);

        double s10 = MS(10, 1, 0);

        double gl = MGL(3, 1, 0);

        System.out.printf("a 10             %,.16f     %,.16f       %,.16f     %,.16f", t10, p10, s10, gl);

        System.out.println();

        t10 = MT(100, 1, 0);

        p10 = MP(100, 1, 0);

        s10 = MS(100, 1, 0);

        gl = MGL(9, 1, 0);

        System.out.printf("b 100             %,.16f     %,.16f       %,.16f     %,.16f", t10, p10, s10, gl);

        System.out.println();

        t10 = MT(1000, 1, 0);

        p10 = MP(1000, 1, 0);

        s10 = MS(1000, 1, 0);

        gl = MGL(27, 1, 0);

        System.out.printf("c 1000            %,.16f     %,.16f       %,.16f     %,.16f", t10, p10, s10, gl);

        System.out.println();

        t10 = MT(10_000, 1, 0);

        p10 = MP(10_000, 1, 0);

        s10 = MS(10_000, 1, 0);

        gl = MGL(81, 1, 0);

        System.out.printf("d 10000           %,.16f     %,.16f       %,.16f     %,.16f", t10, p10, s10, gl);

        System.out.println();

        t10 = MT(100_000, 1, 0);

        p10 = MP(100_000, 1, 0);

        s10 = MS(100_000, 1, 0);

        gl = MGL(243, 1, 0);

        System.out.printf("e 100000          %,.16f     %,.16f       %,.16f     %,.16f", t10, p10, s10, gl);

        System.out.println();

        t10 = MT(1_000_000, 1, 0);

        p10 = MP(1_000_000, 1, 0);

        s10 = MS(1_000_000, 1, 0);

        gl = MGL(729, 1, 0);

        System.out.printf("f 1000000         %,.16f     %,.16f       %,.16f     %,.16f", t10, p10, s10, gl);

        System.out.println();

        t10 = MT(10_000_000, 1, 0);

        p10 = MP(10_000_000, 1, 0);

        s10 = MS(10_000_000, 1, 0);

        gl = MGL(2187, 1, 0);

        System.out.printf("g 10000000        %,.16f     %,.16f       %,.16f     %,.16f", t10, p10, s10, gl);

        System.out.println();

        t10 = MT(100_000_000, 1, 0);

        p10 = MP(100_000_000, 1, 0);

        s10 = MS(100_000_000, 1, 0);

        gl = MGL(6561, 1, 0);

        System.out.printf("h 100000000       %,.16f     %,.16f       %,.16f     %,.16f", t10, p10, s10, gl);

        System.out.println();

        t10 = MT(1_000_000_000, 1, 0);

        s10 = MS(1_000_000_000, 1, 0);

        p10 = MP(1_000_000_000, 1, 0);

        gl = MGL(19683, 1, 0);

        System.out.printf("i 1000000000      %,.16f     %,.16f       %,.16f     %,.16f", t10, p10, s10, gl);

        System.out.println();
    }

    public static double MT(double n, double b, double a) {

        double calka = (zad2(a) + zad2(b)) / 2.0;

        double h = (b - a) / n;

        for (int i = 1; i < n; i++) {

            calka += zad2(a + i * h);

        }


        return calka * h;
    }

    public static double MP(double n, double b, double a) {
        double calka = 0;

        double lew;

        double pra;

        for (double i = 1; i <= n; i++) {

            lew = b - a;
            lew /= 2 * n;

            pra = i;
            pra *= b - a;
            pra /= n;

            calka += zad2(a - lew + pra);
        }

        calka *= (b - a);
        calka /= n;

        return calka;
    }

    public static double MS(double n, double b, double a) {
        double calka = 0;

        double lew;

        double pra;

        for (int i = 1; i <= n; i++) {

            lew = a - (b - a) / (2.0 * n) + i * ((b - a) / n);

            pra = a + i * ((b - a) / n);

            lew = 4 * zad2(lew);

            pra = 2 * zad2(pra);

            calka += pra + lew;
        }

        calka += zad2(a) - zad2(b);

        calka *= (b - a) / (6 * n);

        return calka;
    }

    public static double MGL(double n, double b, double a) {
        double[] x = {-0.906179845938663992797626878299392965125651910762530862873762286,
                -0.538469310105683091036314420700208804967286606905559956202231627,
                0,
                0.538469310105683091036314420700208804967286606905559956202231627,
                0.906179845938663992797626878299392965125651910762530862873762286};
        double[] w = {0.2369268850561890875142640407199173626432600022124140155828278882
                , 0.4786286704993664680412915148356381929122955533431415399727276673
                , 0.5688888888888888888888888888888888888888888888888888888888888888
                , 0.4786286704993664680412915148356381929122955533431415399727276673
                , 0.2369268850561890875142640407199173626432600022124140155828278882};

        double calka = 0, lew, pra, suma = 0;

        for (int i = 1; i <= 5; i++) {
            for (int j = 1; j <= n; j++) {

                lew = (b - a) / (2 * n);

                lew *= (2 * j - 1);

                pra = (b - a) / (2 * n);

                pra *= x[i - 1];

                suma += zad2(a + lew + pra);
            }

            suma *= w[i - 1];

            calka += suma;

            suma = 0;
        }

        calka *= (b - a) / (2 * n);

        return calka;
    }

    public static double MGLPRZE(double N, double b, int a) {
        double x[] = new double[5];
        x[0] = -(1.0 / 3.0) * Math.sqrt(5 + 2 * Math.sqrt(10d / 7d));
        x[1] = -(1.0 / 3.0) * Math.sqrt(5 - 2 * Math.sqrt(10d / 7d));
        x[2] = 0;
        x[3] = (1.0 / 3.0) * Math.sqrt(5 - 2 * Math.sqrt(10d / 7d));
        x[4] = (1.0 / 3.0) * Math.sqrt(5 + 2 * Math.sqrt(10d / 7d));

        double w[] = new double[5];
        w[0] = (322 - 13 * Math.sqrt(70)) / 900;
        w[1] = (322 + 13 * Math.sqrt(70)) / 900;
        w[2] = 128d / 225d;
        w[3] = (322 + 13 * Math.sqrt(70)) / 900;
        w[4] = (322 - 13 * Math.sqrt(70)) / 900;

        double result = 0;
        double calka = 0;
        double temp = (b - a) / N;
        for (int i = 1; i <= 5; i++) {
            result = 0;
            for (int j = 1; j <= N; j++) {
                result += zad2(a + (2 * j - 1) * (temp / 2) + (temp / 2) * x[i - 1]);
            }
            calka += w[i - 1] * result;
        }
        return temp / 2 * calka;
    }


    public static double[] monte() {

        Random los = new Random();

        double y;

        double a = 0, b = 1, suma = 0;

        long n = 10;

        double wynik[] = new double[10];

        for (int i = 0; i < 10; i++) {

            for (long j = 1; j <= n; j++) {
                y = los.nextDouble(1);
                suma += zad2(y);
            }

            suma *= (b - a) / n;

            System.out.println(n + " " + suma);

            wynik[i] = suma;

            n *= 10;
        }


        return wynik;
    }

    public static void monteBlad() {
        double[][] pomiar = new double[10][10];

        double[] srednia = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

        for (int j = 0; j < 10; j++) {
            double z[] = monte();

            for (int i = 0; i < 10; i++) {
                pomiar[i][j] = z[i];
                System.out.print(pomiar[i][j] + " ");
            }
            System.out.println();
        }

        for (int j = 0; j < 10; j++) {
            for (int i = 0; i < 10; i++) {
                srednia[j] += pomiar[j][i];
            }
            srednia[j] /= 10;
            System.out.println(srednia[j]);
        }
    }

    public static void pi() {

        Random los = new Random();

        double m = 0, k = 0, n = 0, nMax = 10;
        double x, y;

        for (int i = 0; i < 7; i++) {

            while (n <= nMax) {
                x = los.nextDouble(1);
                y = los.nextDouble(1);

                if (x * x + y * y <= 1) {
                    k++;
                }

                n++;
            }

            System.out.println(4 * k / n);

            nMax *= 10;

        }
    }

    public static int randomWalker1D(int N) {

        Random los = new Random();

        int pozycja = 0;

        for (int i = 0; i < N; i++) {
            double x = los.nextDouble(1);

            if (x >= 0.5) {
                pozycja++;
            } else if (x < 0.5) {
                pozycja--;
            }
        }

        return pozycja;

    }

    public static void x() {
        int n = 10;
        for (int i = 0; i < 7; i++) {
            System.out.println(randomWalker1D(n));
            n *= 10;
        }
    }

    public static void y() {
        int n = 10;
        double d;
        for (int i = 0; i < 7; i++) {
            double[] x = randomWalker2D(n);
            d = Math.sqrt((x[0] * x[0]) + (x[1] * x[1]));
            System.out.println(d);
            n *= 10;
        }
    }

    public static double[] randomWalker2D(int N) {
        Random los = new Random();

        double x = 0, y = 0;

        for (int i = 0; i < N; i++) {
            double kat = los.nextDouble(1);

            kat *= 2 * Math.PI;

            x += Math.cos(kat);

            y += Math.sin(kat);
        }
        double wynik[] = {x, y};
        return wynik;
    }

    public static void licz() {
        double[] x = {2879.419397, 1881.368332, 2882.912191, 1128.743872,
                3240.884445, 5502.242487, 1904.136855, 3175.941231,
                4287.549899, 1892.034935, 2877.523364

        };
        double suma = 0;
        double srednia = 3.224024465;

        for (int i = 0; i < x.length; i++) {
            suma += Math.pow(x[i] - srednia, 2);
        }

        suma /= 10;

        suma = Math.sqrt(suma);

        System.out.println(suma / 3);

    }


    public static void OH(double tpocz, double tkon, double n) {

        n = 100;

        double h = (tkon - tpocz) / n;

        double x = 0, v = 1, t = tpocz, a = -x;

        for (int i = 0; i < n; i++) {

            a = -x;

            x = x + h * v;

            v = v + h * a;

            t = t + h;

            //System.out.println(x);

            System.out.println(v);

        }

    }

    public static void OHzad2(double tpocz, double tkon, double n) {

        n = 100;

        double h = (tkon - tpocz) / n;

        double x = 0, v = 1, t = tpocz, a = -x;
        for (int j = 0; j < 3; j++) {
            for (int i = 0; i < n; i++) {

                a = -x;

                x = x + h * v;

                v = v + h * a;

                t = t + h;

                System.out.println(x);

                //System.out.println(v);

            }
            n *= 10;
            System.out.println();
            System.out.println();
            System.out.println();
            System.out.println();

        }

    }

    public static void Verlet() {

        double n = 100, tpocz = 0, tkon = 100 * Math.PI;

        double h = (tkon - tpocz) / n;

        double x = 0, v = 1, t = tpocz, a, aNext = -x;
        for (int j = 0; j < 3; j++) {
            for (int i=0;i<n;i++) {

                a = aNext;

                x = x + h * v + Math.pow(h, 2) / 2 * a;

                aNext = -x;

                v = v + h / 2 * (a + aNext);

                t += h;

                //System.out.println(t);

                //System.out.println(x);

                //System.out.println(v);

                System.out.println(aNext);

            }
            n *= 10;
            x = 0;
            v = 1;
            t = tpocz;
            a = -x;
            aNext = -x;
            h = (tkon - tpocz) / n;
            System.out.println();
            System.out.println();
            System.out.println();
            System.out.println();

        }
    }

    public static void lab8() throws InterruptedException {
        int n = 1_000;

        double x = 2.5, v = 1, a, aNext, L = 5, tmin = 0, tmax = 10, h = (tmax - tmin) / n;


        int currentLoc = 150;

        a = 4 * ((12 / Math.pow(x,13) - 6 / Math.pow(x, 7)) - (12 / (Math.pow(L - x, 13)) - 6 / Math.pow(L - x, 7)));

        for(int i=0;i<n;i++){
            x = x + h * v + (h * h) / 2 * a;

            aNext = 4 * ((12 / Math.pow(x, 13) - 6 / Math.pow(x,7)) - (12 / Math.pow(L - x, 13) - 6 / Math.pow(L - x,7)));

            v = v + h / 2 * (a + aNext);

            a = aNext;

        }
    }


    public static void wiz(){

    }
}