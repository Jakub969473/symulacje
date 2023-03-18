import java.util.Scanner;

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

    public static void zadanie1(){

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

    public static double metodaTrapezów(double xp, double xk, double fPocz,
                                        double fKonc){
        /*double dx, calka;
        int n=10000;

        dx = (xk - xp) / (double)n;

        calka = 0;
        for (int i=1; i<n; i++) {
            calka += func(xp + i * dx, pow);
        }
        calka += (func(xp,pow) + func(xk,pow)) / 2;
        calka *= dx;*/

        double calka=xk-xp;

        calka *=(fKonc+fPocz);

        calka /=2;

        return calka;
    }

    public static double metodaProstokątów(double xp, double xk,double fC){

        double calka=xk-xp;

        calka *=fC;

        return calka;
    }

    public static double Simpson(double xp, double xk,double fP,double fK,double fC){

        double calka=(xk-xp)/6;

        double calka2=fK+4*fC+fP;

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

        System.out.println("  N               MT                  MP");

        double t10=MT(10,1,0);

        double p10=MP(10,1,0);

        System.out.println("a 10              "+t10+"  "+p10);

        double t100=MT(100,1,0);

        p10=MP(100,1,0);

        System.out.println("b 100             "+t100+"  "+p10);

        double t1000=MT(1000,1,0);

        p10=MP(1000,1,0);

        System.out.println("c 1000            "+t1000+"  "+p10);

        double t10000=MT(10_000,1,0);

        p10=MP(10_000,1,0);

        System.out.println("d 10000           "+t10000+"  "+p10);

        double t100000=MT(100_000,1,0);

        p10=MP(100_000,1,0);

        System.out.println("e 100000          "+t100000+"  "+p10);

        double t1000000=MT(1_000_000,1,0);

        p10=MP(1_000_000,1,0);

        System.out.println("f 1000000         "+t1000000+" "+p10);

        double t10000000=MT(10_000_000,1,0);

        p10=MP(10_000_000,1,0);

        System.out.println("g 10000000        "+t10000000+" "+p10);

        double t100000000=MT(100_000_000,1,0);

        p10=MP(100_000_000,1,0);

        System.out.println("h 100000000       "+t100000000+"  "+p10);

        double t1000000000=MT(1_000_000_000,1,0);

        p10=MP(1_000_000_000,1,0);

        System.out.println("i 1000000000      "+t1000000000+"  "+p10);
    }

    public static double MT(int n,double b,double a){

        double dx,calka;

        dx = (b - a) / (double)n;

        calka = 0;
        for (int i=1; i<n; i++) {
            calka += zad2(a + i * dx);
        }
        calka += (zad2(a) + zad2(b)) / 2;
        calka *= dx;

        return calka;
    }

    public static double MP(int n,double b,double a){
        double dx,calka=0;

        double lew;

        double pra;

        for (int i=1; i<n; i++) {

            lew=(b-a)/(2*n);

            pra=i*(b-a)/n;

            calka += zad2(a-lew+pra);
        }

        calka*=(b-a)/n;

        return calka;
    }

    public static double MS(){
        double calka=0;

        return calka;
    }

    public static double MGL(){
        double calka=0;

        return calka;
    }
}