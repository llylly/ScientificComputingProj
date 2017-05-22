//
// Created by lly on 21/05/2017.
//

#ifndef SCIENTIFICCOMPUTINGPROJ_PROB3_H
#define SCIENTIFICCOMPUTINGPROJ_PROB3_H

#include "AbstractProb.h"

/**
 * Problem 3: Unit 3 Q6
 * Entrance Function is 'work()'
 * Limyik Li
 */
class Prob3: public AbstractProb {
public:
    int work() override;

private:
    void test(int n);
    double **genHilbert(int n);
    void clean(double **mat, int n, int m);
    int cholesky(double **mat, int n, double ***U, double ***L);
    void printMat(double **mat, int n, int m);
    void solveByChol(double **U, double **L, double *b, int n, double **pX);
};


#endif //SCIENTIFICCOMPUTINGPROJ_PROB3_H
