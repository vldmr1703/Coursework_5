//package functions;
//
//public class TestSystem extends System {
//
//    public TestSystem() {
//        super(2);
//
//    }
//
//    public double[] getX0() {
//        return new double[]{0, 4};
//    }
//
//    public double[] getFunctions(double[] x) {
//        return new double[]{
//                x[0] * 2 + x[1],
//                x[0] * 3 + x[1] * 4
//        };
//    }
//
//    public double[][] getGradient(double[] x) {
//        return new double[][]{
//                new double[]{2, 1},
//                new double[]{3, 4}
//        };
//    }
//
//    public double[] getExactSol(double t) {
//        return new double[]{-Math.exp(t) + Math.exp(5 * t), Math.exp(t) + 3 * Math.exp(5 * t)};
//    }
//
//    public double[][] getM() {
//        return new double[][]{
//                new double[]{1, 0},
//                new double[]{0, 1}
//        };
//    }
//}
