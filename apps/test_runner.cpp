#include <iostream>
#include "symsparseprob.h"

int main()
{
    std::cout << std::endl << "------------------------------------\n\
    SYMMETRIC SPARSE SOLVER TEST\n\
    ------------------------------------\n" << std::endl;

    std::cout << "Test Problem [A][X] = [B]\n\n\
     |  1  5 -3 |       |  8 |\n\
[A]= |  5  2  2 |; [B]= | -2 |\n\
     | -3  2  1 |       |  4 |\n" << std::endl;

    std::cout << "\n------------------------------------\n\nSolver Printout...\n\n" << std::endl;

    unsigned long i = 3;
    SymSparseProb A(i);
    A.setMaxiter(10000);
    A.setToler(1e-15);

    A.setMatrixValue(0, 0, 1., true);
    A.setMatrixValue(1, 1, 2., true);
    A.setMatrixValue(2, 2, 1., true);
    A.setMatrixValue(0, 1, 5., true);
    A.setMatrixValue(0, 2, -3., true);
    A.setMatrixValue(1, 2, 2., true);

    A.B[0] = 8.;
    A.B[1] = -2.;
    A.B[2] = 4.;

    A.JCGSolve();

    A.printOut();

    return 0;
}
