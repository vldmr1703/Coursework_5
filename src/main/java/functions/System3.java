package functions;

public class System3 extends System {
    public System3() {
        super(5);
    }

    @Override
    public double[] getX0() {
        return new double[]{1, 1, 1, 1, 1};
    }

    @Override
    public double[] getFunctions(double[] x) {
        return new double[]{
                2 * x[0] * x[1] * x[2] * x[3],
                -x[0] * x[1] * x[3] * x[3],
                (x[0] * x[1] + x[2] * x[3]) * x[4],
                -x[0] * x[1] * x[1] * x[3] * x[3] * x[4],
                x[0] * x[1] * x[1] - 1
        };
    }

    @Override
    public double[][] getGradient(double[] x) {
        return new double[][]{
                new double[]{
                        2 * x[1] * x[2] * x[3],
                        2 * x[0] * x[2] * x[3],
                        2 * x[0] * x[1] * x[3],
                        2 * x[0] * x[1] * x[2],
                        0
                },
                new double[]{
                        -x[1] * x[3] * x[3],
                        -x[0] * x[3] * x[3],
                        0,
                        -x[0] * x[1] * x[3] * 2,
                        0
                },
                new double[]{
                        x[1] * x[4],
                        x[0] * x[4],
                        x[3] * x[4],
                        x[2] * x[4],
                        x[0] * x[1] + x[2] * x[3]
                },
                new double[]{
                        -x[1] * x[1] * x[3] * x[3] * x[4],
                        -x[0] * 2 * x[1] * x[3] * x[3] * x[4],
                        0,
                        -x[0] * x[1] * x[1] * 2 * x[3] * x[4],
                        -x[0] * x[1] * x[1] * x[3] * x[3]
                },
                new double[]{
                        x[1] * x[1],
                        x[0] * 2 * x[1],
                        0,
                        0,
                        0
                }
        };
    }

    @Override
    public double[] getExactSol(double t) {
        return new double[]{Math.exp(2 * t), Math.exp(-t), Math.exp(2 * t), Math.exp(-t), Math.exp(t)};
    }

    @Override
    public double[][] getM() {
        return new double[][]{
                new double[]{1, 0, 0, 0, 0},
                new double[]{0, 1, 0, 0, 0},
                new double[]{0, 0, 1, 0, 0},
                new double[]{0, 0, 0, 1, 0},
                new double[]{0, 0, 0, 0, 0}
        };
    }

    @Override
    public double getLastDerivative(double[] x) {
        return 0;
    }
}
