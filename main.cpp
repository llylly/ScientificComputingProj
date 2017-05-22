#include <iostream>

#include "Prob1.h"
#include "Prob2.h"
#include "Prob3.h"
#include "Prob4.h"
#include "Prob5.h"
#include "Prob6.h"

#include <iostream>

int main() {
    std::cerr << "Prob1" << std::endl;
    Prob1 *prob1 = new Prob1();
    prob1->work();
    std::cerr << "Prob2" << std::endl;
    Prob2 *prob2 = new Prob2();
    prob2->work();
    std::cerr << "Prob3" << std::endl;
    Prob3 *prob3 = new Prob3();
    prob3->work();
    std::cerr << "Prob4" << std::endl;
    Prob4 *prob4 = new Prob4();
    prob4->work();
    std::cerr << "Prob5" << std::endl;
    Prob5 *prob5 = new Prob5();
    prob5->work();
    std::cerr << "Prob6" << std::endl;
    Prob6 *prob6 = new Prob6();
    prob6->work();
    return 0;
}