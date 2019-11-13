package functions;

public class System4 extends System {

    public System4() {
        super(2);
    }

    @Override
    public double[] getX0() {
        return new double[]{Math.sin(0), -Math.cos(0)};
    }

    @Override
    public double[] getFunctions(double[] x) {
        return new double[]{
        };
    }

    @Override
    public double[] getFunctions(double[] x, double t) {
        return new double[]{
                -x[1], x[0] - Math.sin(t)
        };
    }

    @Override
    public double[][] getGradient(double[] x) {
        return null;
    }

    @Override
    public double[][] getGradient(double[] x, double t) {
        return new double[][]{
                new double[]{
                        0, -1
                },
                new double[]{
                        1 - Math.cos(t), - Math.cos(t)
                }
        };
    }

    @Override
    public double[] getExactSol(double t) {
        return new double[]{Math.sin(t), -Math.cos(t)};
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
        return 0;
    }

    @Override
    public double getLastDerivative(double[] x, double t) {
        return - Math.cos(t);
    }
}
