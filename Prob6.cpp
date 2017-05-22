//
// Created by lly on 22/05/2017.
//

#include <cstdio>
#include <cmath>
#include "Prob6.h"

#define ABS(x) (((x) < 0.0)? (-(x)) : (x))

int Prob6::work() {
    n = 15;
    dataX = new double[15]
            {1, 1.5, 2,
             2.5, 3.0, 3.5,
             4, 4.5, 5,
             5.5, 6, 6.5,
             7, 7.5, 8};
    dataY = new double[15]
            {33.40, 79.50, 122.65,
             159.05, 189.15, 214.15,
             238.65, 252.2, 267.55,
             280.50, 296.65, 301.65,
             310.40, 318.15, 325.15};
    fitting(1);
    fitting(2);
    return 0;
}

void Prob6::fitting(int level) {
    double **samples = new double*[level + 1];
    for (int i=0; i<=level; ++i) {
        samples[i] = getSample(i, dataX, n);
    }
    double **mat = new double*[level + 1];
    double *b = new double[level + 1];
    for (int i=0; i<=level; ++i) {
        mat[i] = new double[level + 1];
        for (int j=0; j<=level; ++j) {
            mat[i][j] = 0.0;
            for (int k=0; k<n; ++k)
                mat[i][j] += samples[i][k] * samples[j][k];
        }
        b[i] = 0;
        for (int k=0; k<n; ++k)
            b[i] += samples[i][k] * dataY[k];
    }
    double *x = solve(mat, b, level + 1);
    printf("Fitting by %d power func:\n", level);
    printf("  y = ");
    for (int i=level; i>=0; --i) {
        printf("(%lf) ", x[i]);
        if (i > 0)
            printf("* x^%d + ", i);
    }
    printf("\n");

    double sd = 0.0;
    for (int i=0; i<n; ++i) {
        double est = 0.0;
        double tmp = 1.0;
        for (int j=0; j<=level; ++j) {
            est += x[j] * tmp;
            tmp *= dataX[i];
        }
        printf("   X: %lf, est: %lf, real: %lf\n", dataX[i], est, dataY[i]);
        sd += (est - dataY[i]) * (est - dataY[i]);
    }
    sd = sqrt(sd / (double)n);
    printf("  RMSE = %lf\n", sd);

    for (int i=0; i<=level; ++i) delete[] samples[i];
    delete[] samples;
    for (int i=0; i<=level; ++i) delete[] mat[i];
    delete[] mat;
    delete[] b;
    delete[] x;
}

double *Prob6::getSample(int level, double *dataX, int n) {
    double *arr = new double[n];
    for (int i=0; i<n; ++i) {
        double tmp = 1;
        for (int j=0; j<level; ++j)
            tmp *= dataX[i];
        arr[i] = tmp;
    }
    return arr;
}

double *Prob6::solve(double **A, double *b, int n) {
    double *x = new double[n];
    for (int i=0; i<n; ++i) {
        // Select pivot
        int p = i;
        for (int j=i+1; j<n; ++j)
            if (ABS(A[j][i]) > ABS(A[p][i]))
                p = j;
        // Exchange pivot line
        if (p != i) {
            double t;
            for (int j = i; j<n; ++j) {
                t = A[i][j], A[i][j] = A[p][j], A[p][j] = t;
            }
            t = b[i], b[i] = b[p], b[p] = t;
        }
        // Flush
        for (int j=i+1; j<n; ++j) {
            double v = A[j][i] / A[i][i];
            for (int k=i; k<n; ++k)
                A[j][k] = A[j][k] - A[i][k] * v;
            b[j] = b[j] - b[i] * v;
        }
    }
    for (int i=n-1; i>=0; --i) {
        x[i] = b[i] / A[i][i];
        for (int j=0; j<i; ++j)
            b[j] -= A[j][i] * x[i];
    }
    return x;
}