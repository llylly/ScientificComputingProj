//
// Created by lly on 21/05/2017.
//

#include <cmath>
#include <cstdio>
#include "Prob3.h"

#define ABS(x) (((x) < 0.0)? (-(x)) : (x))

int Prob3::work() {
    test(8);
    test(10);
    test(12);
    return 0;
}

void Prob3::test(int n) {
    double **mat = genHilbert(n);

    double ***pU = new double**, ***pL = new double**;
    if (cholesky(mat, n, pU, pL)) {
        printf("  Not a positive matrix.\n");
        return;
    }
    double **U = *pU, **L = *pL;
    double *b = new double[n];
    for (int i=0; i<n; ++i) {
        double t = 0.0;
        for (int j=0; j<n; ++j)
            t += mat[i][j];
        b[i] = t;
    }
    double **pX = new double*;
    solveByChol(U, L, b, n, pX);
    double *X = *pX;

    double *r = new double[n];
    for (int i=0; i<n; ++i) {
        r[i] = 0.0;
        for (int j=0; j<n; ++j)
            r[i] += mat[i][j] * X[j];
        r[i] = b[i] - r[i];
    }

    double rMax = 0.0f, deltaxMax = 0.0f;
    for (int i=0; i<n; ++i)
        if (ABS(r[i]) > rMax) rMax = ABS(r[i]);
    for (int i=0; i<n; ++i)
        if (ABS(X[i] - 1.0) > deltaxMax) deltaxMax = ABS(X[i] - 1.0);
    printf("n = %d\n", n);
    printf("  ||r||_inf = %.8e\n", rMax);
    printf("  ||\\delta X||_inf = %.8lf\n", deltaxMax);

    double *bNew = new double[n];
    for (int i=0; i<n; ++i)
        bNew[i] = b[i];
    bNew[0] += 1e-7;
    double **pXnew = new double*;
    solveByChol(U, L, bNew, n, pXnew);
    double *Xnew = *pXnew;

    double *rNew = new double[n];
    for (int i=0; i<n; ++i) {
        rNew[i] = 0.0;
        for (int j=0; j<n; ++j)
            rNew[i] += mat[i][j] * Xnew[j];
        rNew[i] = b[i] - rNew[i];
    }

    double rMaxNew = 0.0f, deltaxMaxNew = 0.0f;
    for (int i=0; i<n; ++i)
        if (ABS(rNew[i]) > rMaxNew) rMaxNew = ABS(rNew[i]);
    for (int i=0; i<n; ++i)
        if (ABS(Xnew[i] - 1.0) > deltaxMaxNew) deltaxMaxNew = ABS(Xnew[i] - 1.0);
    printf(" For 1e-7 turbulence,\n");
    printf("  ||r||_inf = %.8e\n", rMaxNew);
    printf("  ||\\delta X||_inf = %.8lf\n", deltaxMaxNew);

    delete[] X;
    delete[] b;
    clean(U, n, n);
    clean(L, n, n);
    clean(mat, n, n);
    delete[] r;
    delete[] bNew;
    delete[] Xnew;
    delete[] rNew;
    delete pU;
    delete pL;
    delete pX;
    delete pXnew;
}

double **Prob3::genHilbert(int n) {
    double **mat = new double*[n];
    for (int i=0; i<n; ++i)
        mat[i] = new double[n];
    for (int i=0; i<n; ++i)
        for (int j=0; j<n; ++j)
            mat[i][j] = 1.0 / (double)(i+j+1);
    return mat;
}

void Prob3::clean(double **mat, int n, int m) {
    for (int i=0; i<n; ++i)
        delete[] mat[i];
    delete[] mat;
}

int Prob3::cholesky(double **mat, int n, double ***U, double ***L) {
    double **u = new double*[n];
    double **l = new double*[n];
    for (int i=0; i<n; ++i) {
        u[i] = new double[n];
        l[i] = new double[n];
        for (int j=0; j<n; ++j)
            l[i][j] = mat[i][j];
    }
    *U = u, *L = l;
    for (int j=0; j<n; ++j) {
        for (int k=0; k<j; ++k)
            l[j][j] -= l[j][k]*l[j][k];
        if (l[j][j] < 0) return -1;
        l[j][j] = sqrt(l[j][j]);
        for (int i=j+1; i<n; ++i) {
            for (int k = 0; k < j; ++k)
                l[i][j] -= l[i][k] * l[j][k];
            l[i][j] = l[i][j] / l[j][j];
        }
    }
    for (int i=0; i<n; ++i)
        for (int j=i+1; j<n; ++j)
            l[i][j] = 0.0;
    for (int i=0; i<n; ++i)
        for (int j=0; j<n; ++j)
            if (i>j)
                u[i][j] = 0.0;
            else
                u[i][j] = l[j][i];
    return 0;
}

void Prob3::printMat(double **mat, int n, int m) {
    printf("[\n");
    for (int i=0; i<n; ++i) {
        for (int j=0; j<m; ++j)
            if (j < m-1)
                printf("%lf, ", mat[i][j]);
            else
                printf("%lf;\n", mat[i][j]);
    }
    printf("]\n");
}

void Prob3::solveByChol(double **U, double **L, double *b, int n, double **pX) {
    double *y = new double[n];
    double *tb = new double[n];
    double **u = new double*[n];
    double **l = new double*[n];

    for (int i=0; i<n; ++i) {
        u[i] = new double[n];
        l[i] = new double[n];
        for (int j=0; j<n; ++j)
            u[i][j] = U[i][j], l[i][j] = L[i][j];
        tb[i] = b[i];
    }

    for (int i=0; i<n; ++i) {
        y[i] = tb[i] / l[i][i];
        for (int j=i+1; j<n; ++j)
            tb[j] -= y[i] * l[j][i];
    }

    double *X = new double[n];
    for (int i=n-1; i>=0; --i) {
        X[i] = y[i] / u[i][i];
        for (int j=0; j<i; ++j)
            y[j] -= X[i] * u[j][i];
    }

    clean(u, n, n);
    clean(l, n, n);
    delete[] y;
    delete[] tb;

    *pX = X;
}