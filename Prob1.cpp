//
// Created by lly on 21/05/2017.
//

#include "Prob1.h"
#include <cstdio>
#include <cmath>
#include <iostream>

using namespace std;

#define ABS(x) (((x) < 0.0f) ? (-(x)) : (x))

int Prob1::work() {
    float a1 = 0.0f, pre1 = a1, delta1;
    double a2 = 0.0f, pre2 = a2, delta2;

    /* -------------- */
    int i = 0;
    do {
        ++i;
        pre1 = a1, pre2 = a2;
        a1 += 1.0 / i;
        a2 += 1.0 / i;
        delta1 = a1 - pre1;
        delta2 = a2 - pre2;
    } while (delta1 > 0);

    int n = i;
    printf("When using single float, the sum doesn't change when n = %d\n", n);
    printf("Calculated by \"single result - double result\", Error = %lf, Error rate = %lf\n", a1 - a2, (a1 - a2) / a2);

    /* -------------- */
    clock_t start, end;

    a2 = 0.0f;
    long long limit = 100000000;

    start = clock();

    for (long long j=0; j<limit; ++j) {
        a2 += 1.0 / j;
    }

    end = clock();

    double time = (double)(end - start) / (double)CLOCKS_PER_SEC;

    printf("It costs %lf s for %lld terms\n", time, limit);

    /* -------------- */
    double l = 1.0f, r = 1e+20, mid;
    while (r-l > 1.0f) {
        mid = (l+r) / 2.0f;
        float v = log(mid) / log(exp(1.0));
        double epsL = 0.0f, epsR = 1.0f, epsMid;
        while (epsR - epsL > 1e-20) {
            epsMid = (epsL + epsR) / 2.0f;
            float v1 = v + (epsMid);
            if (v == v1)
                epsL = epsMid;
            else
                epsR = epsMid;
        }
        epsMid = (epsL + epsR) / 2.0f;
        if (epsMid * mid < 1.0f) l = mid; else r = mid;
    }
    printf("Estimate unchange N for single = %lld\n", (long long)mid);

    /* -------------- */
    l = 1.0f, r = 1e+20;
    while (r-l > 1.0f) {
        mid = (l+r) / 2.0f;
        double v = log(mid) / log(exp(1.0));
        double epsL = 0.0f, epsR = 1.0f, epsMid;
        while (epsR - epsL > 1e-20) {
            epsMid = (epsL + epsR) / 2.0f;
            double v1 = v + (epsMid);
            if (v == v1)
                epsL = epsMid;
            else
                epsR = epsMid;
        }
        epsMid = (epsL + epsR) / 2.0f;
        if (epsMid * mid < 1.0f) l = mid; else r = mid;
    }
    printf("Estimate unchange N for double = %lld\n", (long long)mid);

    /* -------------- */
    double estN = (double)((long long)mid);
    double estTime = estN * time / limit;
    printf("Therefore, estimate time = %lf s\n", estTime);

    return 0;
}