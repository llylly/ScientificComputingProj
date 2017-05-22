//
// Created by lly on 22/05/2017.
//

#ifndef SCIENTIFICCOMPUTINGPROJ_PROB5_H
#define SCIENTIFICCOMPUTINGPROJ_PROB5_H

#include "AbstractProb.h"

/**
 * Problem 5: Unit 5 Q1
 * Entrance Function is 'work()'
 * Limyik Li
 */
class Prob5: public AbstractProb {
public:
    int work() override;
private:
    void findMaxEigen(double **mat, int n, double *v0, double eps);

    double **mat1, **mat2;
    int dim1, dim2;
    double *mul(double **A, double *b, int n);
    double dot(double *a, double *b, int n);
    double getMaxEle(double *a, int n);
    void clean(double **mat, int n, int m);
};


#endif //SCIENTIFICCOMPUTINGPROJ_PROB5_H
