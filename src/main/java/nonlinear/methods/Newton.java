package nonlinear.methods;

import system.solver.Gauss;

public class Newton {
    protected Function function;
    private static final double EPS = 0.00001;

    protected void print(double[] x) {
        for (double xi : x) {
            System.out.print(xi + " ");
        }
        System.out.println();
    }

    protected boolean comparing(double[] x) {
        double[] F = function.getFunction(x);
        for (int i = 0; i < x.length; i++) {
            if (Math.abs(F[i]) >= EPS)
                return true;
        }
        return false;
    }

    public Newton(Function function) {
        this.function = function;
    }

    public double[] calculate() {
        Gauss system;
        double[] x0 = function.getX0();
        int n = x0.length;
        double[] x1 = new double[n];
        double[] F;
        double[] mas;
        for (int i = 0; i < n; i++) {
            x1[i] = x0[i];
        }
        do {
            for (int i = 0; i < n; i++) {
                x0[i] = x1[i];
            }
            F = function.getFunction(x0);

            system = new Gauss(function.getJacobian(x0), F);
            mas = system.getX();
            for (int i = 0; i < n; i++) {
                x1[i] = x0[i] - mas[i];
            }
            System.out.print("x0 = ");
            print(x0);
            System.out.print("x1 = ");
            print(x1);
        }
        while (comparing(x1));
        return x1;
    }
}

