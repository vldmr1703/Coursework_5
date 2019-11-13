package nonlinear.methods;

import Jama.Matrix;

import functions.Function;
import functions.System;
import system.solver.Gauss;

import static java.lang.System.out;

public class Method2 {

    protected System system;
    protected int n;
    protected static final double EPS_IN = 0.00001;
    protected int m;
    public double h;

    public Method2(System system, double h, int m) {
        this.system = system;
        this.h = h;
        this.m = m;
        this.n = system.getN();
    }

    public double[][] calculate() {
        double t0 = 0;
        int k = 0;
        double[] x0 = getX0();
        double[][] res = new double[m][n];
        res[0] = x0.clone();
        double[] F;
        double[] mas;
        double[] x1;
        do {
            k++;
            double[] previousY = new double[]{0, 0};
            for (int i = 1; i < m; i++) {
                if (i != 1) {
                    previousY = res[i - 2].clone();
                }
                int u = 0;
                x1 = res[i - 1].clone();
                do {
                    u++;
                    x0 = x1.clone();
                    F = getFunctions(x0.clone(), previousY.clone());
                    mas = new Gauss(getGradient(x0.clone()), F.clone()).getX();
//                    mas = new Matrix(getGradient(x0.clone())).inverse().times(new Matrix(F.clone(), 1).transpose()).transpose().getArray()[0];
                    for (int j = 0; j < n; j++) {
                        x1[j] = x0[j] - mas[j];
                    }
                    java.lang.System.out.print("x0 = ");
                    print(x0);
                    java.lang.System.out.print("x1 = ");
                    print(x1);
                } while (u < 2);
                res[i] = x0.clone();
                out.println("u = " + u);
            }
        }
        while (k < 1);
//        while (notEnoughSmall(y));
//        while (mNorm(res[1], res1[1]) > EPS_IN);
        return res;
    }

    public double[] getX0() {
        return new double[]{Math.sqrt(2) / 2, Math.sqrt(2) / 2};
//        return new double[]{0, 4};
    }

    public double[] getFunctions(double[] x, double[] previous) {
        return new double[]{(x[0] - previous[0]) - h * x[1], x[0] * x[0] + x[1] * x[1] - 1};
//        return new double[]{
//                x[0] * 2 + x[1] - (x[0] - previous[0]) / h,
//                x[0] * 3 + x[1] * 4 - (x[1] - previous[1]) / h
//        };
    }

    public double[][] getGradient(double[] x) {
        return new double[][]{
                new double[]{
                        1, -h
                },
                new double[]{
                        2 * x[0], 2 * x[1]
                }
        };
//        return new double[][]{
//                new double[]{2 - 1.0 / h, 1},
//                new double[]{3, 4 - 1.0 / h}
//        };
    }

    public double[] getExactSol(double t) {
        return new double[]{Math.sin(t + Math.PI / 4), Math.cos(t + Math.PI / 4)};
//        return new double[]{-Math.exp(t) + Math.exp(5 * t), Math.exp(t) + 3 * Math.exp(5 * t)};
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
