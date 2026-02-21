
#include <cstdlib>
#include <iostream>
#include "symsparseprob.h"


using namespace std;

// TEST CLASS PER IL SOLUTORE CG

int main()
{
	cout<<endl<<"------------------------------------\n\
    SYMMETRIC SPARSE SOLVER TEST\n\
------------------------------------\n\n";
	cout<<"Test Problem [A][X] = [B]\n\n\
     |  1  5 -3 |       |  8 |\n[A]= |  5  2  2 |; [B]= | -2 |\n     | -3  2  1 |       |  4 |\n\n";
	cout<<"\n------------------------------------\n\nSolver Printout...\n\n";


    unsigned long i = 3;
    SymSparseProb A(i);
    A.setMaxiter(10000);
    A.setMatrixValue(0,0,1.,true);
    A.setMatrixValue(1,1,2.,true);
    A.setMatrixValue(2,2,1.,true);
    A.setMatrixValue(0,1,5.,true);
    A.setMatrixValue(0,2,-3.,true);
    A.setMatrixValue(1,2,2.,true);
	
	A.B[0] = 8.;    
    A.B[1] = -2.;
	A.B[2] = 4.;

    //A.CGSolve();
    A.JCGSolve();
    //A.SSORCGSolve();

    A.printOut();

    return 0;
}
