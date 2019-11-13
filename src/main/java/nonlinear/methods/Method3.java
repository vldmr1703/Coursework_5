package nonlinear.methods;

import Jama.Matrix;

import functions.System;
import system.solver.Gauss;

public class Method3 {

    protected System system;
    protected int n;
    protected static final double EPS_IN = 0.00001;
    protected int m;
    public double h;
    private double a;
    private double b;

    public Method3(System system, double h, int m, double a, double b) {
        this.system = system;
        this.h = h;
        this.m = m;
        this.n = system.getN();
        this.a = a;
        this.b = b;
    }

//    public double[][] calculate() {
////        double[] x = new double[n];
//        double[][] y = new double[m][n];
//        double[] u = new double[n];
//        double[][] A = new double[n][n];
//        double[][] gradient = null;
//        double[][] D = new double[n][n];
//        for (int j = 0; j < n; j++) {
//            D[j][j] = 1;
//        }
//        double[] Fy = new double[n];
//        double[] g = new double[n];
//        double[] b = new double[n];
//        double[] f = new double[n];
//        double lastVal;
//        y[0] = system.getX0();
//        g = system.getFunctions(system.getX0(), 0);
//        for (int j = 0; j < n; j++) {
//            y[1][j] += h * g[j];
//        }
//        g = system.getFunctions(y[1].clone(), 0);
//        for (int i = 1; i < m; i++) {
//            java.lang.System.out.println("i = " + i);
//            while (notEnoughSmall(g)) {
//                u = y[i].clone();
//                lastVal = y[i][n - 1];
////                y[i][n - 1] = g[n - 1];
//                gradient = system.getGradient(y[i-1].clone(), (i-1) * h);
//                double[][] FyD = new double[n][n];
//                for (int j = 0; j < n; j++) {
//                    FyD[j] = gradient[j].clone();
//                    gradient[j] = new double[n];
//                    gradient[j][n - 1] = FyD[j][n - 1];
//                    FyD[j][n - 1] = 0;
//                }
//                FyD = new Matrix(FyD).times(new Matrix(D)).getArray();
////                for (int j = 0; j < n; j++) {
////                    for (int k = 0; k < n - 1; k++) {
////                        gradient[j][k] += FyD[j][k];
////                    }
////                }
//                gradient = new Matrix(gradient).plus(new Matrix(FyD)).getArray();
////                g = system.getFunctions(x.clone());
////                for (int j = 0; j < n - 1; j++) {
////                    g[j] = 0;
////                }
//                b = new Matrix(gradient).inverse().times(new Matrix(g, 1).transpose()).transpose().getArray()[0];
//                y[i][n - 1] = lastVal - b[n - 1];
//                f = system.getFunctions(y[i].clone(), i * h);
//                for (int j = 0; j < n - 1; j++) {
//                    y[i][j] += h * f[j];
//                }
//                g = system.getFunctions(y[i].clone(), i * h);
//                D = new Matrix(D)
//                        .plus(
//                                new Matrix(y[i].clone(), 1).transpose()
//                                        .minus(new Matrix(u.clone(), 1).transpose())
//                                        .minus(new Matrix(D).times(new Matrix(b.clone(), 1).transpose()))
//                                        .times(new Matrix(b.clone(), 1))
//                                        .times(new Matrix(b.clone(), 1)
//                                                .times(new Matrix(b.clone(), 1).transpose()).getArray()[0][0]))
//                        .getArray();
//                print(y[i]);
//            }
//        }
//        return y;
//    }

    public double[][] calculate() {
        double[][] delta = new double[m][n];
        double[][] res = new double[m][n];
        double[][] res1 = new double[m][n];
        res = makeInitialVector(system.getX0(), system.getFunctions(system.getX0(), 0));
        double[][] grad;
        double[] F;
        double t;
        double[] left;
        double[] right;
        double[] diff;
        double[][] M = system.getM();
        double[] FPrev;
        int k = 0;
        do {
            delta = new double[m][n];
//            for (int i = 1; i < m; i++){
//                res1[i] = res[i].clone();
//            }
            for (int i = 1; i < m; i++) {
                t = a + i * h;
                grad = system.getGradient(res[i - 1], t);
                F = system.getFunctions(res[i - 1], t);

                //get last delta
                double num = F[n - 1];
                for (int j = 0; j < n - 1; j++) {
                    num += grad[n - 1][j] * delta[i - 1][j];
                }
                delta[i - 1][n - 1] = -num / grad[n - 1][n - 1];

                //main procedure
//                left = multiplyMatrixVector(grad, delta[i - 1].clone()).clone();
//                right = vectorMinusVector(F.clone(), multiplyMatrixVector(M, getDiff(res[i - 1].clone(), res[i].clone()))).clone();
//                double[] k1 = getF(grad, delta[i - 1], F, M, res[i - 1], res[i]);
//
//                for (int j = 0; j < n - 1; j++) {
//                    delta[i][j] = delta[i - 1][j] + h * (k1[j]);
//                }

                double[] k1 = getF(grad, delta[i - 1], F, M, res[i - 1], res[i]);
                double[] k2 = new double[n];
                double[] k3 = new double[n];
                double[] k4 = new double[n];
                double[] deltaH = new double[n];
                for (int j = 0; j < n; j++) {
                    deltaH[j] = delta[i - 1][j] + h * k1[j] / 2;
                }
                k2 = getF(grad, deltaH, F, M, res[i - 1], res[i]);
                for (int j = 0; j < n; j++) {
                    deltaH[j] = delta[i - 1][j] + h * k2[j] / 2;
                }
                k3 = getF(grad, deltaH, F, M, res[i - 1], res[i]);
                for (int j = 0; j < n; j++) {
                    deltaH[j] = delta[i - 1][j] + h * k3[j];
                }
                k4 = getF(grad, deltaH, F, M, res[i - 1], res[i]);
                for (int j = 0; j < n; j++) {
                    delta[i][j] = delta[i - 1][j] + h * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]) / 6;
                }


            }
            for (int i = 0; i < m; i++) {
                for (int j = 0; j < n; j++) {
                    res[i][j] += delta[i][j];
                }
            }
            print(res);
            k++;
        }
//        while (lNorm(delta) > EPS_IN);
        while (k < 1);
        return res;
    }

    public double[] getF(double[][] grad, double[] delta, double[] F, double[][] M, double[] res, double[] res1) {
        double[] left = multiplyMatrixVector(grad, delta.clone()).clone();
        double[] right = vectorMinusVector(F.clone(), multiplyMatrixVector(M, getDiff(res.clone(), res1.clone()))).clone();
        double[] v = new double[n];
        for (int j = 0; j < n; j++) {
            v[j] = right[j] + left[j];
        }
        return v;
    }

    public double[][] makeInitialVector(double[] x0, double[] F0) {
        double[][] f = new double[m][n];
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                f[i][j] = x0[j] + (a + i * h) * F0[j] + 0.00001 * Math.sin(a+i*h);
            }
//            f[i][0] = Math.sin(a + i * h + Math.PI / 4) + 0.01;
//            f[i][1] = Math.cos(a + i * h + Math.PI / 4) + 0.01;

        }
        return f;
    }

    public double[] getDiff(double[] FPrev, double[] F) {
        double[] diff = new double[n];
        for (int j = 0; j < n; j++) {
            diff[j] = (F[j] - FPrev[j]) / h;
        }
        return diff;
    }

    public double[] getDiff1(double[] res, double[] F) {
        double[] diff = new double[n];
        for (int j = 0; j < n; j++) {
            diff[j] = res[j] + h * F[j];
        }
        return diff;
    }

    public double[] vectorMinusVector(double[] a, double[] b) {
        double[] diff = new double[n];
        for (int j = 0; j < n; j++) {
            diff[j] = a[j] - b[j];
        }
        return diff;
    }


    public double[][] calculate1Index() {
//        double[] x = system.getX0();
        double[][] y = new double[m][n];
        double[] f = new double[n];
        y[0] = system.getX0();
        f = system.getFunctions(y[0].clone(), 0);
        double[][] gradient;
        while (notEnoughSmall(f[n - 1])) {
            y[0][n - 1] -= f[n - 1] / system.getLastDerivative(y[0], 0);
        }
        for (int i = 1; i < m; i++) {
            f = system.getFunctions(y[i - 1].clone(), i * h);
            for (int j = 0; j < n; j++) {
                y[i][j] = y[i - 1][j] + h * f[j];
            }
            f = system.getFunctions(y[i].clone(), i * h);
            int k = 0;
            while (notEnoughSmall(f[n - 1]))
//                while (k<1)

            {
                double div = system.getLastDerivative(y[i].clone(), i * h);
                y[i][n - 1] -= f[n - 1] / div;
                f = system.getFunctions(y[i].clone(), i * h);
                k++;
            }
        }
        print(y);
        return y;
    }

    public double[][] calculate2Index() {
//        double[] x = system.getX0();
        double[][] y = new double[m][n];
        double[] f = new double[n];
        y[0] = system.getX0();
        f = system.getFunctions(y[0].clone(), 0);
        double[][] gradient;
//        while (notEnoughSmall(f[n - 1])) {
//            gradient = system.getGradient(y[0].clone());
//            y[0][n - 1] -= f[n - 1] / system.getLastDerivative(y[0], 0);
//        }
        for (int i = 1; i < m; i++) {
            y[i - 1][0] = Math.sin((i - 1) * h);
            java.lang.System.out.println("i = " + i);
            f = system.getFunctions(y[i - 1].clone(), i * h);
            for (int j = 0; j < n; j++) {
                y[i][j] = y[i - 1][j] + h * f[j];
            }
//            y[i - 1][1] = -y[i][0];

//            y[i][n - 1] = 2;
//            f = system.getFunctions(y[i].clone(), i * h);
//            int k = 0;
//            while (notEnoughSmall(f[n - 1]) && k < 5) {
//                double div = system.getLastDerivative(y[i].clone(), i * h);
//                y[i][n - 1] -= f[n - 1] / div;
//                f = system.getFunctions(y[i].clone(), i * h);
//                k++;
//            }
//            y[i][n - 1] = Math.sqrt(1 - y[i][0]);
//            y[i][0] = Math.sin(i * h);
        }
        return y;
    }

    public double[][] fillInitialMatrix(double[] x0) {
        double[][] functions = new double[m][n];
        double[] f = system.getFunctions(x0);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                functions[i][j] = x0[j]
//                        + i * h * f[j]
                ;
            }
        }
        return functions;
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

    protected double lNorm(double[][] m) {
        double max = 0;
        for (int i = 0; i < m[0].length; i++) {
            double sum = 0;
            for (int j = 0; j < n; j++) {
                sum += Math.abs(m[i][j]);
            }
            if (sum > max) {
                max = sum;
            }
        }
        return max;
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


    protected boolean notEnoughSmall(double[] x) {
//        for (int i = 0; i < x.length; i++) {
        if (Math.abs(x[n - 1]) >= EPS_IN) {
            java.lang.System.out.println(x[n - 1]);
            return true;
//            }
        }
        return false;
    }

    protected boolean notEnoughSmall(double x) {
        if (Math.abs(x) >= EPS_IN) {
            java.lang.System.out.println(x);
            return true;
        }
        return false;
    }

    public double[] multiplyMatrixVector(double[][] A, double[] b) {
        double[] mas = new double[n];
        for (int k = 0; k < n; k++) {
            double s = 0;
            for (int t = 0; t < n; t++) {
                s += A[k][t] * b[t];
            }
            mas[k] = s;
        }
//        return mas;
        return new Matrix(A).times(new Matrix(b, 1).transpose()).transpose().getArray()[0];
    }
}
