package functions;

public class System2 extends System {

    public System2() {
        super(2);
    }

    @Override
    public double[] getX0() {
        return new double[]{
                Math.sqrt(2) / 2,
                Math.sqrt(2) / 2
        };
    }

    @Override
    public double[] getFunctions(double[] x) {
        return new double[]{x[1], x[0] * x[0] + x[1] * x[1] - 1};
    }

    @Override
    public double[] getFunctions(double[] x, double t) {
        return this.getFunctions(x);
    }

    @Override
    public double[][] getGradient(double[] x, double t) {
        return this.getGradient(x);
    }

    @Override
    public double getLastDerivative(double[] x, double t) {
        return this.getLastDerivative(x);
    }

    @Override
    public double[][] getGradient(double[] x) {
        return new double[][]{
                new double[]{
                        0, 1
                },
                new double[]{
                        2 * x[0], 2 * x[1]
                }
        };
    }

    @Override
    public double[] getExactSol(double t) {
        return new double[]{Math.sin(t + Math.PI / 4), Math.cos(t + Math.PI / 4)};
    }

    @Override
    public double[][] getM() {
        return new double[][]{
                new double[]{1, 0},
                new double[]{0, 0}
        };
    }

    @Override
    public double getLastDerivative(double[] x) {
        return 2 * x[1];
    }
}
