//
// Created by lly on 22/05/2017.
//

#include "Prob4.h"

#include <cstdio>

#define ABS(x) (((x) < 0.0)? (-(x)) : (x))

int Prob4::work() {
    double **mat = genHilbert(10);
    double *b = new double[10];
    for (int i=0; i<10; ++i)
        b[i] = 1.0 / (double)(i+1);
    double *X1 = new double[10], *X2 = new double[10];
    for (int i=0; i<10; ++i) X1[i] = X2[i] = 0.0;
    double **pX1 = &X1, **pX2 = &X2;

    int times;

    if ((times = jacobi(mat, b, 10, 1e-4, pX1)) >= 0) {
        printf("Jacobi\n");
        printf(" times = %d\n", times);
        printf(" answer:\n");
        for (int i=0; i<10; ++i)
            printf("  x_%d = %lf\n", i, X1[i]);
        double err = 0.0f;
        for (int i=0; i<10; ++i)
            if (i == 0) {
                if (ABS(X1[i] - 1.0f) > err)
                    err = ABS(X1[i] - 1.0f);
            } else
                if (ABS(X1[i]) > err)
                    err = ABS(X1[i]);
        printf(" ||DeltaX||_inf = %lf\n", err);
    } else
        printf("Jacobi doesn't converge.\n");

    SORTest(mat, b, 10, 1e-4, pX2, 1.25, true);

    for (double oTest = 1e-3; oTest <= 10000000.0; oTest *= 2.0) {
        SORTest(mat, b, 10, 1e-4, pX2, oTest, false);
    }
    for (double oTest = 0.5; oTest <= 1.5; oTest += 0.01) {
        SORTest(mat, b, 10, 1e-4, pX2, oTest, false);
    }

    clean(mat, 10, 10);
    return 0;
}

void Prob4::SORTest(double **mat, double *b, int n, double epsilon, double **pX, double omega, bool printAns) {
    int times;
    for (int i=0; i<n; ++i) (*pX)[i] = 0.0;
    if ((times = SOR(mat, b, n, epsilon, pX, omega)) >= 0) {
        printf("SOR w = %lf\n", omega);
        printf(" times = %d\n", times);
        double err = 0.0f;
        for (int i=0; i<10; ++i)
            if (i == 0) {
                if (ABS((*pX)[i] - 1.0f) > err)
                    err = ABS((*pX)[i] - 1.0f);
            } else
            if (ABS((*pX)[i]) > err)
                err = ABS((*pX)[i]);
        printf(" ||DeltaX||_inf = %lf\n", err);
        if (printAns) {
            printf(" answer:\n");
            for (int i = 0; i < 10; ++i)
                printf("  x_%d = %lf\n", i, (*pX)[i]);
        }
    } else
        printf("SOR doesn't converge when w = %lf.\n", omega);
}

double **Prob4::genHilbert(int n) {
    double **mat = new double*[n];
    for (int i=0; i<n; ++i)
        mat[i] = new double[n];
    for (int i=0; i<n; ++i)
        for (int j=0; j<n; ++j)
            mat[i][j] = 1.0 / (double)(i+j+1);
    return mat;
}

void Prob4::clean(double **mat, int n, int m) {
    for (int i=0; i<n; ++i)
        delete[] mat[i];
    delete[] mat;
}

int Prob4::jacobi(double **mat, double *b, int n, double epsilon, double **pX) {
    double *newX = new double[n];
    double *X = *pX;
    double nowEps;
    int times = 0;
    do {
        ++times;
        nowEps = 0.0;
        for (int i=0; i<n; ++i) {
            newX[i] = b[i];
            for (int j=0; j<n; ++j)
                if (j != i)
                    newX[i] -= mat[i][j] * X[j];
            newX[i] /= mat[i][i];
            if (ABS(newX[i] - X[i]) > nowEps)
                nowEps = ABS(newX[i] - X[i]);
        }
        for (int i=0; i<n; ++i)
            X[i] = newX[i];
    } while ((nowEps >= epsilon) && (times <= 100));
    delete[] newX;
    if (nowEps >= epsilon)
        return -1;
    else
        return times;
}

int Prob4::SOR(double **mat, double *b, int n, double epsilon, double **pX, double omega) {
    double *newX = new double[n];
    double *X = *pX;
    double nowEps;
    int times = 0;
    do {
        ++times;
        nowEps = 0.0;
        for (int i=0; i<n; ++i) {
            newX[i] = b[i];
            for (int j=0; j<i; ++j)
                newX[i] -= mat[i][j] * newX[j];
            for (int j=i+1; j<n; ++j)
                newX[i] -= mat[i][j] * X[j];
            newX[i] /= mat[i][i];
            newX[i] = (1.0 - omega) * X[i] + omega * newX[i];
        }
        for (int i=0; i<n; ++i)
            if (ABS(newX[i] - X[i]) > nowEps)
                nowEps = ABS(newX[i] - X[i]);
        for (int i=0; i<n; ++i)
            X[i] = newX[i];
    } while ((nowEps >= epsilon) && (times <= 10000));
    delete[] newX;
    if (nowEps >= epsilon)
        return -1;
    else return times;
}