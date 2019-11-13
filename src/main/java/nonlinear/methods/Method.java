package nonlinear.methods;

import Jama.Matrix;

import functions.Function;
import functions.System;

public class Method {

    protected System system;
    protected int n;
    protected static final double EPS_IN = 0.00001;
    protected int m;
    public double h;

    public Method(System system, double h, int m) {
        this.system = system;
        this.h = h;
        this.m = m;
        this.n = system.getN();
    }

    public double[][] calculate() {
        double t0 = 0;
        int k = 0;
        double[] x0 = system.getX0();
        double[][] res = null;
        double[][] res1 = fillInitialMatrix(x0);
        double[][] y = new double[m][n];
        double[][] gradient = null;
        double[] f = null;
        do {
            k++;
            double[][] M = system.getM();
//            double[] y0;
            y[0] = new double[n];
            res = copyMatrix(res1);
            for (int i = 1; i < m; i++) {
//                y0 = y[i - 1].clone();
//                y0 = new double[]{0, 0};
                gradient = system.getGradient(res[i - 1].clone());
                f = system.getFunctions(res[i - 1].clone());
//                if(i!=1) {
                for (int j = n - 1; j >= 0; j--) {
                    if (M[j][j] == 0) {
                        double sum = 0;
                        for (int e = 0; e < n; e++) {
                            if (e != j) {
                                sum += gradient[j][e] * y[i - 1][e];
                            }
                        }
                        if(gradient[j][j]==0){
                            gradient[j][j]=1e-9;
                        }
                        y[i - 1][j] = (f[j] - sum) / gradient[j][j];
//                            y[j] = y[0];
//                        res1[i - 1][j] = res[i - 1][j] + y[i - 1][j];
                    }
                }
//                }
                double[] gradF = new Matrix(gradient).times(new Matrix(y[i - 1].clone(), 1).transpose()).transpose().getArray()[0];
                for (int j = n - 1; j >= 0; j--) {
                    if (M[j][j] != 0) {
                        y[i][j] = y[i - 1][j] + h * (gradF[j] -
                                (res[i][j] - res[i - 1][j]) / h
//                                (res[i - 1][j] + h * f[j])
                                + f[j]);
//                        res1[i][j] = res[i][j] + y[i][j];
                    } else {
//                        y[j] = y0[j];
//                        double sum = 0;
//                        for (int e = 0; e < n; e++) {
//                            if (e != j) {
//                                sum += (gradient[j][e]) * y[e];
//                            }
//                        }
//                        y[j] = (f[j] - sum) / gradient[j][j];
////                        res[i - 1][j] = sum;
//                        res1[i - 1][j] = res[i - 1][j] + y[j];
                    }
                }
//                for (int j = n - 1; j >= 0; j--){
//                    if (M[j][j] == 0){
//                        y[j] = y0[j];
//                    }
//                }
            }
            for (int j = n - 1; j >= 0; j--) {
                if (M[j][j] == 0) {
                    double sum = 0;
                    for (int e = 0; e < n; e++) {
                        if (e != j) {
                            sum += (gradient[j][e]) * y[m - 2][e];
                        }
                    }
                    sum = (f[j] - sum )/ gradient[j][j];
                    y[m - 1][j] = sum;
//                    res1[m - 1][j] = res[m - 1][j] + sum;
                }
            }
            for (int i = 0; i < m; i++) {
                for (int j = 0; j < n; j++) {
                    res1[i][j] = res[i][j] + y[i][j];
                }
            }
            print(res1);
//            java.lang.System.out.println(notEnoughSmall(y));
        }
        while (k < 4);
//        while (notEnoughSmall(y));
//        while (mNorm(res[1], res1[1]) > EPS_IN);
        return res;
    }

    public double[][] fillInitialMatrix(double[] x0) {
        double[][] functions = new double[m][n];
        double[] f = system.getFunctions(x0);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                functions[i][j] = x0[j]
                        + i * h * f[j]
                ;
            }
        }
        return functions;
    }

    public Function[] getF(Function[] y, Function[] f) {
        Function[] res = new Function[n];
        for (int i = 0; i < n; i++) {
            res[i] = y[i].add(f[i].mult(h)).add(f[i].mult(-1));
        }
        return res;
    }

    public Function[] getF(Function[] y, Function[] y0, Function[] f) {
        Function[] res = new Function[n];
        for (int i = 0; i < n; i++) {
            res[i] = (y[i].add(y0[i].mult(-1))).mult(1.0 / h).add(f[i].mult(-1));
        }
        return res;
    }


    public void printF(double[] x) {
        for (int i = 1; i < x.length; i++) {
            java.lang.System.out.print("& $" + x[i] + "$ ");
            java.lang.System.out.println("\\");
            java.lang.System.out.println("\\hline");
        }
    }

    public void print(double[] x) {
        for (double xi : x) {
            java.lang.System.out.print(xi + " ");
        }
        java.lang.System.out.println();
    }

    protected double mNorm(double[] x0, double[] x1) {
        double s = 0;
        for (int i = 0; i < x0.length; i++) {
            s += Math.pow(x0[i] - x1[i], 2);
        }
        return Math.sqrt(s);
    }

    protected double error(double[] x0, double[] x1) {
        double max = Math.abs(x0[0] - x1[0]);
        for (int i = 1; i < n; i++) {
            double abs = Math.abs(x0[i] - x1[i]);
            if (abs > max) {
                max = abs;
            }
        }
        return max;
    }

    protected double[][] divDifferences(double[] x0, double[] x1, double t) {
        int n = x0.length;
        double[][] jac = new double[n][n];
        double[] F0;
        double[] F1;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                double[] mas0 = new double[n];
                for (int k = 0; k < j + 1; k++) {
                    mas0[k] = x0[k];
                }
                for (int k = j + 1; k < n; k++) {
                    mas0[k] = x1[k];
                }
                double[] mas1 = new double[n];
                for (int k = 0; k < j; k++) {
                    mas1[k] = x0[k];
                }
                for (int k = j; k < n; k++) {
                    mas1[k] = x1[k];
                }
//                F0 = system.getFunctions(mas0, t);
//                F1 = system.getFunctions(mas1, t);
                double r = (x0[j] - x1[j]);
                if (r == 0) {
                    r = 1e-8;
                }
//                jac[i][j] = (F0[i] - F1[i]) / (r);
            }
        }
        return jac;
    }

    public void print(double[][] A) {
        for (int i = 0; i < A.length; i++) {
            java.lang.System.out.println();
            for (int j = 0; j < A[0].length; j++) {
                java.lang.System.out.print(A[i][j] + " ");
            }
            java.lang.System.out.println();
        }
        java.lang.System.out.println("---------------------------------");
    }

    protected double[][] copyMatrix(double[][] a) {
        double[][] m = new double[a.length][];
        for (int i = 0; i < a.length; i++) {
            m[i] = a[i].clone();
        }
        return m;
    }

    public double norm(double[][] a, double[][] b) {
        double max = 0;
        for (int i = 0; i < n; i++) {
            double sum = 0;
            for (int j = 0; j < m; j++) {
                sum += Math.abs(a[j][i] - b[j][i]);
            }
            if (sum > max) {
                max = sum;
            }
        }
        return max;
    }

    protected boolean comparing(double[] x) {
        for (int i = 0; i < x.length; i++) {
            if (Math.abs(x[i]) >= EPS_IN) {
                java.lang.System.out.println(x[i]);
                return true;
            }
        }
        return false;
    }
}
