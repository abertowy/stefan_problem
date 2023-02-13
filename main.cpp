#include "stefanProblem.h"

#include <QApplication>
#include <QtWidgets/QMainWindow>
#include <QtCharts/QChartView>
#include <QtCharts/QLineSeries>

QT_CHARTS_USE_NAMESPACE

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    StefanProblemSolver solver;
    solver.calculate(21);

    QLineSeries *series_21 = new QLineSeries();
    for (int i = 0; i < solver.PTB.size(); ++i) {
        series_21->append(solver.ts[i], solver.PTB[i]);
    }
    series_21->setName("21 pts ("+ QString::number(solver.solveTime) + " ms)");

    solver.clean();
    solver.calculate(51);

    QLineSeries *series_51 = new QLineSeries();
    for (int i = 0; i < solver.PTB.size(); ++i) {
        series_51->append(solver.ts[i], solver.PTB[i]);
    }
    series_51->setName("51 pts ("+ QString::number(solver.solveTime) + " ms)");

    solver.clean();
    solver.calculate(101);

    QLineSeries *series_101 = new QLineSeries();
    for (int i = 0; i < solver.PTB.size(); ++i) {
        series_101->append(solver.ts[i], solver.PTB[i]);
    }
    series_101->setName("101 pts ("+ QString::number(solver.solveTime) + " ms)");

    solver.clean();
    solver.calculate(201);

    QLineSeries *series_201 = new QLineSeries();
    for (int i = 0; i < solver.PTB.size(); ++i) {
        series_201->append(solver.ts[i], solver.PTB[i]);
    }
    series_201->setName("201 pts ("+ QString::number(solver.solveTime) + " ms)");

    QChart *chart = new QChart();
    chart->legend()->setVisible(true);
    chart->legend()->setAlignment(Qt::AlignBottom);
    chart->addSeries(series_21);
    chart->addSeries(series_51);
    chart->addSeries(series_101);
    chart->addSeries(series_201);
    chart->createDefaultAxes();
    chart->axisX()->setTitleText("ts");
    chart->axisY()->setTitleText("PTB");
    chart->setTitle("StefanProblem in 21 val");

    QChartView *chartView = new QChartView(chart);
    chartView->setRenderHint(QPainter::Antialiasing);

    QMainWindow w;
    w.setCentralWidget(chartView);
    w.resize(1200, 600);
    w.show();
    return a.exec();
}
