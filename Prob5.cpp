//
// Created by lly on 22/05/2017.
//

#include "Prob5.h"
#include <cstdio>

#define ABS(x) (((x) < 0.0)? (-(x)) : (x))

int Prob5::work() {
    dim1 = 3;
    mat1 = new double*[3];
    for (int i=0; i<3; ++i) {
        mat1[i] = new double[3];
        if (i==0) {
            mat1[i][0] = 5, mat1[i][1] = -4, mat1[i][2] = 1;
        } else
        if (i==1) {
            mat1[i][0] = -4, mat1[i][1] = 6, mat1[i][2] = -4;
        } else
        if (i==2) {
            mat1[i][0] = 1, mat1[i][1] = -4, mat1[i][2] = 7;
        }
    }
    dim2 = 4;
    mat2 = new double*[4];
    for (int i=0; i<4; ++i) {
        mat2[i] = new double[4];
        if (i==0) {
            mat2[i][0] = 25, mat2[i][1] = -41, mat2[i][2] = 10, mat2[i][3] = -6;
        } else
        if (i==1) {
            mat2[i][0] = -41, mat2[i][1] = 68, mat2[i][2] = -17, mat2[i][3] = 10;
        } else
        if (i==2) {
            mat2[i][0] = 10, mat2[i][1] = -17, mat2[i][2] = 5, mat2[i][3] = -3;
        } else
        if (i==3) {
            mat2[i][0] = -6, mat2[i][1] = 10, mat2[i][2] = -3, mat2[i][3] = 2;
        }
    }

    double* v0_1 = new double[dim1];
    for (int i=0; i<dim1; ++i)
        v0_1[i] = 1.0;

    double* v0_2 = new double[dim2];
    for (int i=0; i<dim2; ++i)
        v0_2[i] = 1.0;

    findMaxEigen(mat1, dim1, v0_1, 1e-5);
    findMaxEigen(mat2, dim2, v0_2, 1e-5);

    clean(mat1, dim1, dim1);
    clean(mat2, dim2, dim2);
    delete[] v0_1;
    delete[] v0_2;
    return 0;
}

void Prob5::findMaxEigen(double **mat, int n, double *v0, double eps) {
    double *u = new double[n];
    for (int i=0; i<n; ++i) u[i] = v0[i];
    double t = getMaxEle(u, n);
    for (int i=0; i<n; ++i) u[i] /= t;
    double ans = 0.0, delta;
    int i = 0;
    do {
        ++i;
        double *v = mul(mat, u, n);
        double nAns = getMaxEle(v, n);
        // Rayleigh Quotient Optimization
        double rAns = dot(u, v, n) / dot(u, u, n);
        delta = rAns - ans;
        ans = rAns;
        delete[] u;
        u = v;
        for (int i=0; i<n; ++i)
            u[i] /= nAns;
        printf("  i = %d, max(v_k) = %lf, lambda = %lf\n", i, nAns, ans);
        printf("    u = [");
        for (int i=0; i<n; ++i) {
            printf(" %lf", u[i]);
            if (i < n-1) printf(","); else printf(" ]^T\n");
        }
    } while (ABS(delta) >= eps);
    printf("Max Eigenvector = %lf\n", ans);
    delete[] u;
}

double *Prob5::mul(double **A, double *b, int n) {
    double *ans = new double[n];
    for (int i=0; i<n; ++i) {
        ans[i] = 0.0f;
        for (int j=0; j<n; ++j)
            ans[i] += A[i][j] * b[j];
    }
    return ans;
}

double Prob5::dot(double *a, double *b, int n) {
    double ans = 0.0;
    for (int i=0; i<n; ++i)
        ans += a[i] * b[i];
    return ans;
}

double Prob5::getMaxEle(double *a, int n) {
    double ans = 0.0;
    for (int i=0; i<n; ++i)
        if (ABS(a[i]) > ans)
            ans = ABS(a[i]);
    return ans;
}

void Prob5::clean(double **mat, int n, int m) {
    for (int i=0; i<n; ++i)
        delete[] mat[i];
    delete[] mat;
}