//package euler;
//
//import functions.Function;
//import functions.System;
//
//public class Euler {
//    //initial values
//    private Function[] x0;
//    private Function[] y0;
//    //step
//    private double h;
//    //count of equations in the system
//    private int n;
//
//    public Euler(Function[] x0, Function[] y0, double h) {
//        this.x0 = x0.clone();
//        this.y0 = y0.clone();
//        this.h = h;
//        this.n = x0.length;
//    }
//
//    /**
//     * Get solution of the system like:
//     * y1' = f(x,y1,y2);
//     * y2' = f(x,y1,y2);
//     *
//     * @param system
//     * @param m      - count of time fragments
//     * @return
//     */
//    public Function[][] getDerivatives(System system, int m) {
//        Function[][] y = new Function[m][n];
//        y[0] = y0.clone();
//        for (int i = 1; i < m; i++) {
//            Function[] mas = system.getFunctions(y[i - 1].clone());
//            double[] yTemp = new double[n];
//            Function[] functions = new Function[n];
//            double t = (i - 1) * h;
//            for (int j = 0; j < n; j++) {
//                yTemp[j] = y[i - 1][j].add(mas[j].mult(h)).getValue(t);
//                int finalJ = j;
//                functions[j] = p -> yTemp[finalJ];
//            }
//            Function[] mas1 = system.getFunctions(functions);
//            for (int j = 0; j < n; j++) {
//                y[i][j] = y[i - 1][j].add(mas[j].add(mas1[j]).mult(h / 2));
//            }
//        }
//        return y;
//    }
//
//    public Function[] getDerivatives1(System system, int m, double t) {
//        Function[] y;
//        y = y0.clone();
//        Function[] mas = system.getFunctions(y.clone());
//        double[] yTemp = new double[n];
//        Function[] functions = new Function[n];
//        for (int j = 0; j < n; j++) {
//            yTemp[j] = y[j].add(mas[j].mult(h)).getValue(t);
//            int finalJ = j;
//            functions[j] = p -> yTemp[finalJ];
//        }
//        Function[] mas1 = system.getFunctions(functions);
//        for (int j = 0; j < n; j++) {
//            y[j] = y[j].add(mas[j].add(mas1[j]).mult(h / 2));
//        }
//        return y;
//    }
//}
