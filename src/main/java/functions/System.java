package functions;

public abstract class System {

    protected int n;

    public System(int n) {
        this.n = n;
    }

    public int getN() {
        return n;
    }

//    public abstract double[] getX0();

    public abstract double[] getX0();

    public abstract double[] getFunctions(double[] x);

    public double[] getFunctions(double[] x, double t) {
        return null;
    }

    ;

    public abstract double[][] getGradient(double[] x);

    public double[][] getGradient(double[] x, double t) {
        return null;
    }

    ;


    public abstract double[] getExactSol(double t);

    /**
     * x'(t) = ...
     * y'(t) = ...
     * x^2(t) + y^2(t) + z^2(t) = 1
     * |
     * |   1 0 0
     * M = 0 1 0
     * 0 0 0
     */
    public abstract double[][] getM();

    public abstract double getLastDerivative(double[] x);

    public double getLastDerivative(double[] x, double t) {
        return 0;
    }

    ;
}
