#ifndef LINEAR_SOLVER_COMPLEX_H_INCLUDED
#define LINEAR_SOLVER_COMPLEX_H_INCLUDED

#include "sparse_matrix_complex.h"
#include <complex>

using namespace std;

class ComplexLinearSolver
{
    public:
        unsigned long number_of_unknowns;
        complex<double> *solution;
        complex<double> *right_hand_side;
        ComplexSparseMatrix matrix;

        ComplexLinearSolver(unsigned long number_of_unknowns, unsigned int max_iterations = 50000,
                          double tolerance = 1e-15);
        ~ComplexLinearSolver();

        void setTolerance(double tolerance);
        void setMaxIterations(unsigned int max_iterations);
        void setMatrixValue(unsigned long row_index, unsigned long column_index,
                          complex<double> value, bool replace);
        void printResult();
        unsigned int getIterations();

        complex<double> dotProduct(complex<double> *vector_one, complex<double> *vector_two);
        void solveWithCOCG();
        void solveWithCOCGAndSSOR(double relaxation_factor = 1.0);

    private:
        void applySSORPreconditioner(complex<double> *input_vector, complex<double> *output_vector, double omega);

        unsigned int maximum_iterations;
        double solver_tolerance;
        unsigned int iterations;
        double residual;
        bool out_of_iterations;
};

#endif
