
/******************************************************************************

                               LINEAR SOLVER (REAL)
                               EMANT 2013

******************************************************************************/

#ifndef LINEAR_SOLVER_REAL_H_INCLUDED
#define LINEAR_SOLVER_REAL_H_INCLUDED

#include "sparse_matrix.h"

using namespace std;

class RealLinearSolver
{
    public:
        unsigned long number_of_unknowns;
        double *solution;
        double *right_hand_side;
        Sparse matrix;

        RealLinearSolver(unsigned long number_of_unknowns, unsigned int max_iterations = 50000,
                      double tolerance = 1e-15, unsigned int iteration_update = 2000);
        ~RealLinearSolver();

        void setTolerance(double tolerance);
        void setMaxIterations(unsigned int max_iterations);
        void setMatrixValue(unsigned long row_index, unsigned long column_index,
                          double value, bool replace);
        void setIterationUpdate(unsigned int iteration_update);
        void printResult();
        unsigned int getIterations();

        double dotProduct(double *vector_one, double *vector_two);

        void solveWithCG();
        void solveWithJCG();
        void solveWithSSOR(double relaxation_factor = 1.0);

    private:
        void applySSORPreconditioner(double *input_vector, double *output_vector, double omega);

        unsigned int maximum_iterations;
        double solver_tolerance;
        unsigned int iteration_update_frequency;
        unsigned int iterations;
        double residual;
        bool out_of_iterations;
};

inline double RealLinearSolver::dotProduct(double *vector_one, double *vector_two)
{
    double dot_result = 0;
    for (unsigned long index = 0; index < number_of_unknowns; index++)
        dot_result += vector_one[index] * vector_two[index];
    return dot_result;
};

inline void RealLinearSolver::applySSORPreconditioner(double *input_vector, double *output_vector, double omega)
{
    unsigned long index;
    double scale = 1.0 / (omega * (2.0 - omega));
    Entry *element;

    // Forward sweep: solve (D + omega*L) * y = r
    // L is strictly lower triangular. For j < i, A[i][j] is stored in row j at column i
    for (index = 0; index < number_of_unknowns; index++)
    {
        double sum = input_vector[index];
        for (unsigned long j = 0; j < index; j++)
        {
            // Search for column 'index' in row j's linked list
            element = matrix.diag[j].first;
            while (element && element->col < index)
                element = element->next;
            if (element && element->col == index)
                sum -= omega * element->val * output_vector[j];
        }
        output_vector[index] = sum / matrix.diag[index].val;
    }

    // Backward sweep: solve (D + omega*L^T) * z = y
    // L^T is upper triangular. For j > i, A[j][i] = A[i][j] stored in row i
    for (index = number_of_unknowns; index > 0; index--)
    {
        double sum = output_vector[index-1];
        element = matrix.diag[index-1].first;
        while (element)
        {
            sum -= omega * element->val * output_vector[element->col];
            element = element->next;
        }
        output_vector[index-1] = sum / matrix.diag[index-1].val;
    }

    // Final scaling: z = z / (omega * (2 - omega))
    for (index = 0; index < number_of_unknowns; index++)
        output_vector[index] *= scale;
}

#endif
