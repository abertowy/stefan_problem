#include "stefanProblem.h"
#include <algorithm>
#include <cmath>
#include <chrono>
#include <iostream>
#include <QDebug>

//Вычисление определителя квадратной матрицы a[n][n]
/*double determinant(std::vector<std::vector<double>> &a, unsigned int n) {
    unsigned int m = n;
    if (m == 0) return 0;
    if (m == 1) return a[0][0];
    if (m == 2) return (a[0][0] * a[1][1] - a[1][0] * a[0][1]);
    bool sign = false; // смена знака определителя. по умолчанию - нет
    double det = 1; // определитель
    double tmp;
    unsigned int x, y;

    // цикл по всей главной диагонали
    for (unsigned int i = 0; i < n; i++) {
        // выносим элемент a[i][i] за определитель
        det *= a[i][i];
        tmp = a[i][i];
        for (x = i; x < m; x++) {
            a[i][x] = a[i][x] / tmp;
        }
        // таким образом a[i][i] теперь равен 1
        // зануляем все элементы стоящие под (i, i)-ым,
        // при помощи вычитания с опр. коеффициентом
        for (y = i + 1; y < n; y++) {
            tmp = a[y][i];
            for (x = i; x < m; x++)
            a[y][x] -= (a[i][x] * tmp);
        }
    }
    if (sign) return det * (-1);
    return det;
};*/

std::vector<std::vector<double>> matrixTransformation(std::vector<std::vector<double>> a) {
    for (size_t i = 0; i < a.size(); ++i) {
        for (size_t j = 0; j < a.size(); ++j) {
            if (i == j)
                continue;
            std::swap(a[i][j], a[j][i]);
        }
    }
    return a;
}

std::vector<double> solveHauss(std::vector<std::vector<double>> a, std::vector<double> &b, const size_t &n) {
    std::vector<double> res;
    for(size_t i = 0; i < n;  ++i) {
        a[i].push_back(b[i]);
    }
    for (size_t rowIdx = 0; rowIdx < n; ++rowIdx) {
        if (!a[rowIdx][rowIdx]) {
            for(size_t rowForSwap = 0; rowForSwap < n; ++rowForSwap) {
                if (rowForSwap == rowIdx)
                    continue;
                if (a[rowForSwap][rowIdx]) {
                    a[rowIdx].swap(a[rowForSwap]);
                    break;
                }
            }
            return res;
        }
        if (a[rowIdx][rowIdx] != 1) {
            for (long long divColumnIdx = n; divColumnIdx >= 0; --divColumnIdx) {
                a[rowIdx][divColumnIdx] = a[rowIdx][divColumnIdx]/a[rowIdx][rowIdx];
            }
        }
        for (size_t incrRow = rowIdx+1; incrRow < n; ++incrRow) {
            for (long long incrColumn = n; incrColumn >= 0; --incrColumn) {
                a[incrRow][incrColumn] = a[incrRow][incrColumn] - a[incrRow][rowIdx]*a[rowIdx][incrColumn];
            }
        }
    }
    for(long long revRowIdx = n-1; revRowIdx >= 0; --revRowIdx) {
        for(long long incrRow = revRowIdx-1; incrRow >= 0; --incrRow) {
            for (long long decrColumn = n; decrColumn >= revRowIdx; --decrColumn)
                a[incrRow][decrColumn] = a[incrRow][decrColumn] - a[incrRow][revRowIdx]*a[revRowIdx][decrColumn];
        }
    }
    res.reserve(n);
    for (size_t resRowIdx = 0; resRowIdx < n; ++resRowIdx)
        res.push_back(a[resRowIdx][n]);
    return res;
}

std::vector<double> multiplyMatrixVector(std::vector<std::vector<double>> &a, std::vector<double> &b, const size_t &n) {
    std::vector<double> res;
    res.reserve(a.size());
    double res_row;

    for(size_t i = 0; i < a.size(); ++i) {
        res_row = 0.0;
        for (size_t j = 0; j < n; ++j) {
            res_row += a[i][j]*b[j];
        }
        res.push_back(res_row);
    }
    return res;
}

/*double determinantGauss(std::vector<std::vector<double>> &a, unsigned int n) {
    //std::vector<std::vector<double>> tA (matrixTransformation(a));
    std::vector<std::vector<double>> re(n, std::vector<double>(n));
    for (int i = 0; i < n ++i)
        re[i][i] = 1;
    for (int columnIdx = 0; columnIdx < n; ++columnIdx) {
        for (int rowIdx = columnIdx; rowIdx < n; ++rowIdx) {
            if (a[rowIdx][columnIdx]) {
                if (rowIdx != columnIdx) {
                    a[columnIdx].swap(a[rowIdx]);
                    for(int j = columnIdx; j < n; ++j) {
                        a[columnIdx][j] = a[columnIdx][j]/a[columnIdx][columnIdx];
                        re[columnIdx][j] = re[columnIdx][j]/a[columnIdx][columnIdx];
                    }
                    for (int decrRow = columnIdx+1; decrRow < n; ++decrRow) {
                        re[decrRow][columnIdx] = -a[decrRow][columnIdx];
                        for (int decrColumn = columnIdx; decrColumn < n; ++decrColumn)
                            a[decrRow][decrColumn] = a[decrRow][decrColumn] - a[decrRow][columnIdx]*a[columnIdx][decrColumn];

                    }
                }
            }
            break;
        }
    }

    for(int revRowIdx = n-1; revRowIdx >= 0; --revRowIdx) {
        for(int incrRow = revRowIdx-1; incrRow >= 0; --incrRow) {
            re[incrRow][revRowIdx] = -a[incrRow][revRowIdx];
            for (int decrColumn = n-1; decrColumn >= revRowIdx; --decrColumn)
                a[incrRow][decrColumn] = a[incrRow][decrColumn] - a[incrRow][revRowIdx]*a[revRowIdx][decrColumn];
        }
    }

}*/

void StefanProblemSolver::clean() {
    PTB.clear();
    PTB.push_back(0.0);
    phi.clear();
    ts.clear();
    ts.push_back(0.0);
    solveTime = 0;
}

void StefanProblemSolver::calculate(unsigned int val) {
    double l = 2.0;
    double dx=l/(val-1);
    std::vector <double> xs, T;
    double it = 0.0;
    for (size_t i = 0; i < val; ++i){
        xs.push_back(it);
        phi.push_back(it - PTB[0]);                  // initialize the signed distance function based of the starting PTB
        T.push_back(phi.back() > 0 ? Tm_ : T_);  // initialize the temperature function to T_ on the left and Tm everywhere that is liquid
        it = it + dx;
    }
    double dt = 0.01/365.0;
    double t = t0_;

    // Write matrices for calculating the heat flux, distance gradient and diffusive terms

    double VN = ks_*dt/(pow(dx,2));
    double VN_reinit = k*kd/(pow(dx,2));
    std::vector<std::vector<double>> heat_flux(val, std::vector<double>(val)), center(val, std::vector<double>(val)), re(val, std::vector<double>(val));
    heat_flux[0][0] = -(Ks_/(dx));
    re[0][0] = 1+2*VN_reinit;
    for(size_t i = 1; i < val; ++i) {
        heat_flux[i][i] = -(Ks_/(dx));
        heat_flux[i-1][i] = Ks_/(dx);

        center[i][i] = 1+2*VN;
        center[i-1][i] = -VN;
        center[i][i-1] = -VN;

        re[i][i] = 1+2*VN_reinit;
        re[i-1][i] = -VN_reinit;
        re[i][i-1] = -VN_reinit;
    }
    center[0][0] = 1;
    center[0][1] = 0;
    center[val-1][val-2] = 0;
    center[val-1][val-1] = 1;

    re[0][1] = -2.0*VN_reinit;
    re[val-1][val-2] = -2.0*VN_reinit;

    std::chrono::time_point<std::chrono::steady_clock> start = std::chrono::steady_clock::now();

    std::vector <double> grad_phi(phi.size(), 0);
    while (t < 0.1) {
        std::vector<double> heatFluxVec (multiplyMatrixVector(heat_flux, T, heat_flux.size()));
        for (size_t i = 0; i < phi.size(); ++i) {
            if (i == 0)
                grad_phi[i] = (phi[i+1] - phi[i])/dx;
            else if (i == phi.size() -1)
                grad_phi[i] = (phi[i] - phi[i-1])/dx;
            else
                grad_phi[i] = (phi[i+1] - phi[i-1])/(2*dx);
        }

        for (size_t i = 0; i < phi.size(); ++i) {
            double speed = 1.0/(rho_*Lf_)*heatFluxVec[i]*(grad_phi[i] < 0.0 ? -1.0 : 1.0);
            phi[i] -= dt*speed*std::fabs(grad_phi[i]);
        }

        reinitPhi(re, dx);
        T = solveHauss(center, T, center.size());

        std::vector<double> reducedXs;
        reducedXs.reserve(xs.size());
        for (size_t i = 0; i < T.size(); ++i) {
            if (phi[i] > 0.0 && T[i] < 0.0)
                T[i] = 0.0;
            if(T[i] >= -0.01)
                reducedXs.push_back(xs[i]);
        }
        t += dt;

        ts.push_back(t);
        PTB.push_back(*std::min_element(reducedXs.begin(), reducedXs.end()));
    }

    std::chrono::time_point<std::chrono::steady_clock> end = std::chrono::steady_clock::now();
    std::chrono::milliseconds diff = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    solveTime = diff.count();
}

void StefanProblemSolver::reinitPhi(std::vector<std::vector<double>> &re_val, double dx) {
    size_t size = phi.size();
    std::vector <double> se, d_phi(size, 1);
    for (auto val: phi) {
        se.push_back(val < 0 ? -1 : 1);
    }

    double count = 0.0;
    auto it = std::max_element(d_phi.begin(), d_phi.end(), [](int a, int b) { return std::fabs(a) < std::fabs(b); });
    double d_phi_max = it != d_phi.end() ? *it : 0;

    while (d_phi_max > 0.01) {
        for (size_t i = 0; i < phi.size(); ++i) {
            if (i == 0)
                d_phi[i] = k*se[i]*(1.0 - std::fabs((phi[i+1] - phi[i])/dx));
            else if (i == phi.size() -1)
                d_phi[i] = k*se[i]*(1.0 - std::fabs((phi[i] - phi[i-1])/dx));
            else
                d_phi[i] = k*se[i]*(1.0 - std::fabs((phi[i+1] - phi[i-1])/(2*dx)));
        }
        for (size_t i = 0; i < size; ++i)
            phi[i] += d_phi[i];

        phi[0] -= k*kd*2/dx;
        phi[size-1] -= -k*kd*2/dx;
        phi = solveHauss(re_val,phi, re_val.size());
        ++count;

        auto it = std::max_element(d_phi.begin(), d_phi.end(), [](int a, int b) { return std::fabs(a) < std::fabs(b); });
        d_phi_max = it != d_phi.end() ? *it : 0;
    }
}
