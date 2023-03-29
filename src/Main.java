import java.util.Scanner;
import java.util.function.Function;

public class Main {
    public static void main(String[] args) {

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

        zadanie2();

        /*przykład1();

        przykład2();

        przykład3();

        przykład4();*/

        //test();

    }

    public static void test() {

        double a = 0.0;

        double b = 1.0;

        int n = 10;

        for (int i = 0; i < 9; i++) {

            double result1 = trapezoidalMethod( a, b, n);

            System.out.println(result1);
            n=n*10;
        }

/// policzenia wyniku całkowania dla każdej z funkcji testowych

        double result = trapezoidalMethod(a,b,n );

    }

    public static double trapezoidalMethod(double a,double b,int n){// Eugodt poprzeszkolu

        double h=(b-a)/n;

        double sum = (zad2(a) + zad2(b))/2.0; // suna wartości funkcji na krancach przedziału

        for (int i = 1; i < n; i++) {
            double x = a + i*h; // punkt środkowy 1-tegu pooprzedziału

            sum+= zad2(x); // dodanie wartotel funkcji w punkcie środkowye de suny

        }

        return h*sum;
        // wynik satuwaele przybliżony, autody Trapezow


}


    public static void xxx(){
       double [] x={ 2.436,
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

       double wynik=0;

       for(int i=0;i<x.length;i++){
           wynik+=Math.pow(x[i]-2.44,2);
       }

       System.out.format("%f",wynik);


    }


    public static void zzz(){
        double [] x={ 2.433,
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

        double wynik=0;

        for(int i=0;i<x.length;i++){
            wynik+=Math.pow(x[i]-2.45,2);
        }

        System.out.format("%f",wynik);


    }
    public static void Float(){
        long millisActualTime = System.currentTimeMillis();

        float a=2;

        float b=3;

        float c=6;

        float d=7;

        float x;

        float y;

        float z;

        float t;

        float u;

        double time=1*Math.pow(10,9);

        for(int i=0;i<time;i++){
            x=a/c;

            y=a/b;

            z=d/c;

            t=(b/d)*c;

            u=(b*c)/d;

            a++;

            b++;

            c++;

            d++;

        }


        long executionTime = System.currentTimeMillis() - millisActualTime;

        System.out.println((double)executionTime/1000);
    }

    public static void Double(){
        long millisActualTime = System.currentTimeMillis();

        double a=2;

        double b=3;

        double c=6;

        double d=7;

        double x;

        double y;

        double z;

        double t;

        double u;

        double time=1*Math.pow(10,9);

        for(int i=0;i<time;i++){
            x=a/c;

            y=a/b;

            z=d/c;

            t=(b/d)*c;

            u=(b*c)/d;

            a++;

            b++;

            c++;

            d++;

        }


        long executionTime = System.currentTimeMillis() - millisActualTime;

        System.out.println((double)executionTime/1000);
    }

    public static void bisekcja(double a,double b){

        //x^4-15x^2-15x+5;

        double eps=0.00000000000001;

        double c=0;

        double funA,funB,funC;

        funA=funkcja(a);

        funB=funkcja(b);

        if(funA*funB>=0){
            throw new RuntimeException("funkcja nie ma miejsc zerowych");
        }

        while(Math.abs(b-a)>eps){

            c=(a+b)/2;

            funC=funkcja(c);

            if(funA*funC<0){
                b=c;
            }else{
                a=c;
                funA=funC;
            }

        }

        System.out.print("Miejsce zerowe to ");

        System.out.format("%.16f",c);

    }

    public static double funkcja(double wartość){

        double v = Math.pow(wartość, 4) - 15 * Math.pow(wartość, 2) - 15 * wartość + 5;
        return v;
    }

    public static void przykład1(){
        System.out.println(metodaTrapezów(0,10,0));

        System.out.println(metodaTrapezów(-1,1,1));

        System.out.println(metodaTrapezów(-10,10,2));

        System.out.println(metodaTrapezów(-10,0,3));

        System.out.println(metodaTrapezów(0,1,4));

    }

    public static void przykład2(){
        System.out.println(metodaProstokątów(0,10,0));

        System.out.println(metodaProstokątów(-1,1,1));

        System.out.println(metodaProstokątów(-10,10,2));

        System.out.println(metodaProstokątów(-10,0,3));

        System.out.println(metodaProstokątów(0,1,4));

    }

    public static void przykład3(){
        System.out.println(Simpson(0,10,0));

        System.out.println(Simpson(-1,1,1));

        System.out.println(Simpson(-10,10,2));

        System.out.println(Simpson(-10,0,3));

        System.out.println(Simpson(0,1,4));

    }
    public static void przykład4(){

        System.out.println("wynik całek    wynik metody Gaus");

        double aT = Gaus(0, 10, 0);

        System.out.println("90              "+aT);

        double bT = Gaus(-1, 1,1);

        System.out.println("-2              "+bT);

        double cT = Gaus(-10, 10,2);

        System.out.println("686.666         "+cT);

        double dT = Gaus(-10, 0,3);

        System.out.println("-10 100         "+dT);

        double eT = Gaus(0, 1,4);

        System.out.printf("0               %,.32f",eT);
        System.out.println();

        double fT = Gaus(0, 1,5);

        System.out.printf("0               %,.32f",fT);
        System.out.println();

        double gT = Gaus(0, 1,6);

        System.out.printf("0               %,.17f",gT);
        System.out.println();

    }

    public static double metodaTrapezów(double xp, double xk, int pow){
        double calka=xk-xp;

        calka *=(func(xp,pow)+func(xk,pow));

        calka /=2;

        return calka;
    }

    public static double metodaProstokątów(double xp, double xk,int pow){

        double calka=xk-xp;

        double x=xp+xk;

        calka *=func(x/2,pow);

        return calka;
    }

    public static double Simpson(double xp, double xk,int pow){

        double calka=(xk-xp)/6;

        double x=xp+xk;

        double calka2=func(xp,pow)+4*func(x/2,pow)+func(xk,pow);

        calka*=calka2;

        return calka;
    }

    public static double Gaus(double xp, double xk,int pow) {
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

        double calka=0;

        double temp;

        double temp2;

        for (int i = 0; i < 5; i++) {

            temp=(xk-xp)/2;

            temp*=x[i];

            temp2=(xk+xp)/2;

            temp=func(temp+temp2,pow);
            
            temp*=w[i];

            calka+=temp;
        }

        calka*=((xk-xp)/2);

        return calka;
    }
    public static double func(double x,int pow){

        if(pow==0){
            return fA(x);
        }else if(pow==1){
            return fB(x);
        }else if(pow==2){
            return fC(x);
        }else if(pow==3){
            return fD(x);
        }else if(pow==4){
            return fE(x);
        }else if(pow==5){
            return fF(x);
        }else if(pow==6){
            return fG(x);
        }

        return 0;
    }

    public static double fA(double x){
        return 9;
    }
    public static double fB(double x){
        return (3*x)-1;
    }
    public static double fC(double x){
        return Math.pow(x,2)+1;
    }
    public static double fD(double x){
        return  4*Math.pow(x,3)+2*x;
    }
    public static double fE(double x){
        return 5*Math.pow(x,4)-2*x;
    }
    public static double fF(double x){
        return 10*Math.pow(x,9)-8*Math.pow(x,7);
    }
    public static double fG(double x){
        return 15*Math.pow(x,14)-12*Math.pow(x,11);
    }

    public static double zad2(double x){return 354*Math.pow(x,117);}

    private static void zadanie2() {

        System.out.println("  N               MT                     MP                       MS                     MGL");

        double t10=MT(10,1,0);

        double p10=MP(10,1,0);

        double s10=MS(10,1,0);

        double gl=MGL(10,1,0);

        System.out.printf("a 10             %,.16f     %,.16f       %,.16f     %,.16f",t10,p10,s10,gl);

        System.out.println();

        t10=MT(100,1,0);

        p10=MP(100,1,0);

        s10=MS(100,1,0);

        gl=MGL(100,1,0);

        System.out.printf("b 100             %,.16f     %,.16f       %,.16f     %,.16f",t10,p10,s10,gl);

        System.out.println();

        t10=MT(1000,1,0);

        p10=MP(1000,1,0);

        s10=MS(1000,1,0);

        gl=MGL(1000,1,0);

        System.out.printf("c 1000            %,.16f     %,.16f       %,.16f     %,.16f",t10,p10,s10,gl);

        System.out.println();

        t10=MT(10_000,1,0);

        p10=MP(10_000,1,0);

        s10=MS(10_000,1,0);

        gl=MGL(10_000,1,0);

        System.out.printf("d 10000           %,.16f     %,.16f       %,.16f     %,.16f",t10,p10,s10,gl);

        System.out.println();

        t10=MT(100_000,1,0);

        p10=MP(100_000,1,0);

        s10=MS(100_000,1,0);

        gl=MGL(100_000,1,0);

        System.out.printf("e 100000          %,.16f     %,.16f       %,.16f     %,.16f",t10,p10,s10,gl);

        System.out.println();

        t10=MT(1_000_000,1,0);

        p10=MP(1_000_000,1,0);

        s10=MS(1_000_000,1,0);

        gl=MGL(1_000_000,1,0);

        System.out.printf("f 1000000         %,.16f     %,.16f       %,.16f     %,.16f",t10,p10,s10,gl);

        System.out.println();

        t10=MT(10_000_000,1,0);

        p10=MP(10_000_000,1,0);

        s10=MS(10_000_000,1,0);

        gl=MGL(10_000_000,1,0);

        System.out.printf("g 10000000        %,.16f     %,.16f       %,.16f     %,.16f",t10,p10,s10,gl);

        System.out.println();

        t10=MT(100_000_000,1,0);

        p10=MP(100_000_000,1,0);

        s10=MS(100_000_000,1,0);

        gl=MGL(100_000_000,1,0);

        System.out.printf("h 100000000       %,.16f     %,.16f       %,.16f     %,.16f",t10,p10,s10,gl);

        System.out.println();

        t10=MT(1_000_000_000,1,0);

        s10=MS(1_000_000_000,1,0);

        p10=MP(1_000_000_000,1,0);

        gl=MGL(1_000_000_000,1,0);

        System.out.printf("i 1000000000      %,.16f     %,.16f       %,.16f     %,.16f",t10,p10,s10,gl);

        System.out.println();
    }

    public static double MT(double n,double b,double a){

        double calka=(zad2(a)+zad2(b))/2.0;

        double h=(b-a)/n;

        for (int i=1; i<n; i++) {

            calka+=zad2(a+i*h);

        }


        return calka*h;
    }

    public static double MP(double n,double b,double a){
        double calka=0;

        double lew;

        double pra;

        for (double i=1; i<=n; i++) {

            lew=b-a;
            lew/=2*n;

            pra=i;
            pra*=b-a;
            pra/=n;

            calka += zad2(a-lew+pra);
        }

        calka*=(b-a);
        calka/=n;

        return calka;
    }

    public static double MS(double n,double b,double a){
        double calka=0;

        double lew;

        double pra;

        for (int i=1; i<=n; i++) {

            lew=a-(b-a)/(2.0*n)+i*((b-a)/n);

            pra=a+i*((b-a)/n);

            lew=4*zad2(lew);

            pra=2*zad2(pra);

            calka += pra+lew;
        }

        calka+=zad2(a)-zad2(b);

        calka*=(b-a)/(6*n);

        return calka;
    }

    public static double MGL(double n,double b,double a){
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

        double calka=0,lew,pra,suma=0;

        for(int i=1;i<=5;i++) {
            for (int j = 1; j <= n; j++) {

                lew = (b - a) / (2 * n);

                lew *= (2 * j - 1);

                pra = (b - a) / (2 * n);

                pra *= x[i-1];

                suma += zad2(a + lew + pra);
            }

            suma*=w[i-1];

            calka+=suma;

            suma=0;
        }

        calka*=(b-a)/(2*n);

        return calka;
    }

    public static double MGLPRZE(double N, double b, int a){
            double x[] = new double[5];
            x[0] = -(1.0/3.0) * Math.sqrt(5 + 2 * Math.sqrt(10d / 7d));
            x[1] = -(1.0/3.0) * Math.sqrt(5 - 2 * Math.sqrt(10d / 7d));
            x[2] = 0;
            x[3] = (1.0/3.0) * Math.sqrt(5 - 2 * Math.sqrt(10d / 7d));
            x[4] = (1.0/3.0) * Math.sqrt(5 + 2 * Math.sqrt(10d / 7d));

            double w[] = new double[5];
            w[0] = (322 - 13 * Math.sqrt(70)) / 900;
            w[1] = (322 + 13 * Math.sqrt(70)) / 900;
            w[2] = 128d / 225d;
            w[3] = (322 + 13 * Math.sqrt(70)) / 900;
            w[4] = (322 - 13 * Math.sqrt(70)) / 900;

            double result = 0;
            double calka=0;
            double temp = (b-a)/N;
            for (int i = 1; i <= 5; i++) {
                result=0;
                for (int j = 1; j <= N; j++) {
                    result += zad2(a + (2*j - 1)* (temp/2) + (temp/2)*x[i-1]);
                }
                calka += w[i-1]*result;
            }
            return temp/2 * calka;
        }
}