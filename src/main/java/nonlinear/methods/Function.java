package nonlinear.methods;

public abstract class Function {

    protected int n;

    public Function(int n) {
        this.n = n;
    }

    public int getN() {
        return n;
    }

    public abstract double[] getX0();

    public abstract double[] getFunction(double[] x);

    public abstract double[][] getJacobian(double[] x0);
}