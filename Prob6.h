//
// Created by lly on 22/05/2017.
//

#ifndef SCIENTIFICCOMPUTINGPROJ_PROB6_H
#define SCIENTIFICCOMPUTINGPROJ_PROB6_H

#include "AbstractProb.h"

/**
 * Problem 6: Unit 6 Q3
 * Entrance Function is 'work()'
 * Limyik Li
 */
class Prob6: public AbstractProb {
public:
    int work() override;
private:
    void fitting(int level, double *dataX, double *dataY, bool lnFitting);

    double *dataX, *dataY;
    int n;
    double *getSample(int level, double *dataX, int n);
    double *solve(double **A, double *b, int n);
};


#endif //SCIENTIFICCOMPUTINGPROJ_PROB6_H
