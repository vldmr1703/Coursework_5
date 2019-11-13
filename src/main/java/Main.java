import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.NumberTickUnit;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.ui.ApplicationFrame;
import org.jfree.data.category.DefaultCategoryDataset;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import functions.System2;
import functions.System3;
import functions.System4;
import nonlinear.methods.Method3;

public class Main {
    public static double f(double x, double y) {
        return x * x + y * y - 1;
    }

    public static double diff(double x, double y) {
        return 2 * y;
    }

    public static void main(String[] args) {
        int n = 2;
        int m = 500;
        double a = 0;
        double b = 6.28;
        double h = (b - a) / m;
//        double t = 1.000741541384761;
//        double t1 = 0.014727737953702025;int k=0;
//        while (Math.abs(f(t, t1)) > 0.0001) {
//            t1 -= f(t, t1) / diff(t, t1);
//            k++;
//        }
//        java.lang.System.out.println(k);
        System2 testSystem = new System2();
        Method3 method = new Method3(testSystem, h, m, a, b);
//        new Matrix(new double[]{1, 2, 3}, 1)
//                .times(new Matrix(new double[]{1, 2, 3}, 1).transpose());
        double[][] system = method.calculate();
//        for (int i = 0; i < m; i++) {
//            for (int j = 0; j < n; j++) {
////                System.out.print(system[i][j].getValue(i * h) + " ");
////                System.out.print(testSystem.getExactSol(j * 0.1)[i] + " ");
//            }
//            for (int j = 0; j < n; j++) {
//                java.lang.System.out.print(testSystem.getExactSol(i * h)[j] + " ");
////                System.out.print(testSystem.getExactSol(j * 0.1)[i] + " ");
//            }
//            java.lang.System.out.println();
//        }

        XYSeriesCollection xySeriesCollection = new XYSeriesCollection();
        XYSeries[] exactSeries = new XYSeries[n];
        for (int j = 0; j < n; j++) {
            exactSeries[j] = new XYSeries("x" + (j + 1));
            for (int i = 0; i < m; i++) {
                exactSeries[j].add(a + i * h, testSystem.getExactSol(i * h)[j]);
            }
        }
        for (int i = 0; i < n; i++) {
            xySeriesCollection.addSeries(exactSeries[i]);
        }

        XYSeries[] approxSeries = new XYSeries[n];
        for (int j = 0; j < n; j++) {
            approxSeries[j] = new XYSeries("approxX" + (j + 1));
            for (int i = 0; i < m; i++) {
                approxSeries[j].add(a + i * h, system[i][j]);
            }
        }
        for (int i = 0; i < n; i++) {
            xySeriesCollection.addSeries(approxSeries[i]);
        }


        XYDataset dataset = xySeriesCollection;
        DrawGraph chart = new DrawGraph(
                "",
                "", dataset, a, b);

        chart.pack();
        chart.setVisible(true);
    }

    static class DrawGraph extends ApplicationFrame {

        public DrawGraph(String applicationTitle, String chartTitle, XYDataset dataset, double a, double b) {
            super(applicationTitle);
            JFreeChart lineChart = ChartFactory.createXYLineChart(
                    chartTitle,
                    "Time", "Function value",
                    dataset,
                    PlotOrientation.VERTICAL,
                    true, true, false);
            XYPlot xyPlot = lineChart.getXYPlot();
            NumberAxis domain = (NumberAxis) xyPlot.getDomainAxis();
            domain.setRange(a, b);
            domain.setTickUnit(new NumberTickUnit((b - a) / 10));

            NumberAxis range = (NumberAxis) xyPlot.getRangeAxis();
            range.setRange(-1.5, 10);
            range.setTickUnit(new NumberTickUnit(0.5));
            ChartPanel chartPanel = new ChartPanel(lineChart);
            chartPanel.setPreferredSize(new java.awt.Dimension(700, 500));
            setContentPane(chartPanel);
        }

        private DefaultCategoryDataset createDataset() {
            int n = 10;
            DefaultCategoryDataset dataset = new DefaultCategoryDataset();
            double fault = 0;
            return dataset;
        }
    }
}