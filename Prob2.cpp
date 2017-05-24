//
// Created by lly on 21/05/2017.
//

#include "Prob2.h"

#include <cstdio>

#define ABS(x) (((x) < 0.0)? (-(x)) : (x))

int Prob2::work() {
    double x1 = newtonSolve(&Prob2::f1, &Prob2::df1, 0.6, 1.0, 1e-6);
    printf("Final answer for (1): %.8lf\n", x1);
    double x2 = newtonSolve(&Prob2::f2, &Prob2::df2, 1.2, 1.0, 1e-6);
    printf("Final answer for (2): %.8lf\n", x2);
}

double Prob2::newtonSolve(double (*f)(double), double (*df)(double), double x0, double lambda0, double epsilon) {
    double x = x0, preX = x0, nexX = x0;
    double lambda = lambda0;
    while ((ABS(f(x)) > epsilon) || (ABS(x - preX) > epsilon)) {
        double s = f(x) / df(x);
        nexX = x - s;
        while (ABS(f(nexX)) >= ABS(f(x))) {
            nexX = x - lambda * s;
            lambda /= 2.0;
        }
        preX = x, x = nexX;
        printf("lambda = %.8lf, x = %.8lf, x - pre_x = %.8lf, f(x) = %.8lf\n", lambda, x, x - preX, f(x));
    }
    return x;
}

double Prob2::f1(double x) {
    return x*x*x - x - 1.0;
}

double Prob2::df1(double x) {
    return 3.0*x*x - 1.0;
}

double Prob2::f2(double x) {
    return -x*x*x + 5.0*x;
}

double Prob2::df2(double x) {
    return -3.0*x*x + 5.0;
}

