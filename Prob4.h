//
// Created by lly on 22/05/2017.
//

#ifndef SCIENTIFICCOMPUTINGPROJ_PROB4_H
#define SCIENTIFICCOMPUTINGPROJ_PROB4_H

#include "AbstractProb.h"


/**
 * Problem 4: Unit 4 Q1
 * Entrance Function is 'work()'
 * Limyik Li
 */
class Prob4: public AbstractProb {
public:
    int work() override;

private:
    void SORTest(double **mat, double *b, int n, double epsilon, double **pX, double omega, bool print);

    double **genHilbert(int n);
    void clean(double **mat, int n, int m);
    int jacobi(double **mat, double *b, int n, double epsilon, double **pX);
    int SOR(double **mat, double *b, int n, double epsilon, double **pX, double omega);
};


#endif //SCIENTIFICCOMPUTINGPROJ_PROB4_H
