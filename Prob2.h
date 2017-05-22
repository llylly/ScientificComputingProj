//
// Created by lly on 21/05/2017.
//

#ifndef SCIENTIFICCOMPUTINGPROJ_PROB2_H
#define SCIENTIFICCOMPUTINGPROJ_PROB2_H

#include "AbstractProb.h"

/**
 * Problem 2: Unit 2 Q2
 * Entrance Function is 'work()'
 * Limyik Li
 */
class Prob2: public AbstractProb {
public:
    int work() override;

private:
    double newtonSolve(double (*f)(double), double (*df)(double), double x0, double lambda0, double epsilon);
    static double f1(double x);
    static double df1(double x);
    static double f2(double x);
    static double df2(double x);
};


#endif //SCIENTIFICCOMPUTINGPROJ_PROB2_H
