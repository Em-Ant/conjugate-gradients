
/******************************************************************************

                               LINEAR SOLVER (REAL)
                               EMANT 2013

******************************************************************************/

#include <cmath>
#include "linear_solver_real.h"
#include <iostream>

RealLinearSolver::RealLinearSolver(unsigned long n,
                unsigned int max_iterations, double tolerance, unsigned int iteration_update) : matrix(n)
{
    number_of_unknowns = n;
    solution = new double[number_of_unknowns]();
    right_hand_side = new double[number_of_unknowns]();

    solver_tolerance = tolerance;
    maximum_iterations = max_iterations;
    iteration_update_frequency = iteration_update;
    iterations = 0;
    residual = 0.;
    out_of_iterations = false;
}

RealLinearSolver::~RealLinearSolver()
{
    delete[] solution;
    delete[] right_hand_side;
}

void RealLinearSolver::setTolerance(double tolerance)
{
    solver_tolerance = tolerance;
}

void RealLinearSolver::setMaxIterations(unsigned int max_iterations)
{
    maximum_iterations = max_iterations;
}

void RealLinearSolver::setIterationUpdate(unsigned int iteration_update)
{
    iteration_update_frequency = iteration_update;
}

void RealLinearSolver::setMatrixValue(unsigned long row_index, unsigned long column_index,
                                   double value, bool replace)
{
    matrix.setValue(row_index, column_index, value, replace);
}

unsigned int RealLinearSolver::getIterations()
{
    return iterations;
}

void RealLinearSolver::printResult()
{
    std::cout << " - Linear Solver (Real) Data : -" << std::endl << std::endl;
    std::cout << "[Matrix]";
    matrix.printOut();
    std::cout << "[Right Hand Side]        [Solution]" << std::endl << std::endl;
    for(unsigned long index = 0; index < number_of_unknowns; index++)
        std::cout << "[" << right_hand_side[index] << "]        [" << solution[index] << "]" << std::endl;
    std::cout << std::endl << std::endl;
}

void RealLinearSolver::solveWithCG()
{
    unsigned long row_index, iteration_count;
    double *matrix_times_direction;
    double *search_direction;
    double *residual_vector;
    double alpha, beta;
    double residual_norm_old, residual_norm_new = 0.0;

    residual_vector = new double[number_of_unknowns]();
    search_direction = new double[number_of_unknowns]();
    matrix_times_direction = new double[number_of_unknowns]();
    double right_hand_side_norm = 0.;
    out_of_iterations = false;

    for(row_index = 0; row_index < number_of_unknowns; row_index++)
        right_hand_side_norm += right_hand_side[row_index] * right_hand_side[row_index];
    right_hand_side_norm = std::sqrt(right_hand_side_norm);

    matrix.gaxpy(solution, residual_vector);
    for(row_index = 0; row_index < number_of_unknowns; row_index++)
    {
        residual_vector[row_index] = right_hand_side[row_index] - residual_vector[row_index];
        search_direction[row_index] = residual_vector[row_index];
    }

    residual_norm_old = dotProduct(residual_vector, residual_vector);

    iteration_count = 0;
    do
    {
        residual_norm_new = 0.;
        matrix_times_direction = new double[number_of_unknowns]();
        matrix.gaxpy(search_direction, matrix_times_direction);

        double denominator = dotProduct(search_direction, matrix_times_direction);

        if (std::abs(denominator) < 1e-30)
        {
            std::cout << "CG: denominator zero, stopping" << std::endl;
            break;
        }

        alpha = residual_norm_old / denominator;

        for(row_index = 0; row_index < number_of_unknowns; row_index++)
        {
            solution[row_index] += alpha * search_direction[row_index];
            residual_vector[row_index] -= alpha * matrix_times_direction[row_index];
        }

        residual_norm_new = dotProduct(residual_vector, residual_vector);
        double current_residual = std::sqrt(std::abs(residual_norm_new));

        if (current_residual < solver_tolerance * right_hand_side_norm)
            break;

        beta = residual_norm_new / residual_norm_old;

        for(row_index = 0; row_index < number_of_unknowns; row_index++)
            search_direction[row_index] = residual_vector[row_index] + beta * search_direction[row_index];

        residual_norm_old = residual_norm_new;

        iteration_count++;
        if (iteration_count > maximum_iterations)
        {
            out_of_iterations = true;
            break;
        }
        delete[] matrix_times_direction;
    } while (true);

    delete[] residual_vector;
    delete[] search_direction;
    delete[] matrix_times_direction;

    iterations = iteration_count;
    residual = std::sqrt(residual_norm_new);
}

void RealLinearSolver::solveWithJCG()
{
    unsigned long row_index, iteration_count;
    double *matrix_times_direction;
    double *search_direction;
    double *residual_vector;
    double *preconditioned_residual;
    double alpha, beta;
    double residual_norm_old, residual_norm_new = 0.0;
    double *inverse_diagonal;

    residual_vector = new double[number_of_unknowns]();
    search_direction = new double[number_of_unknowns]();
    matrix_times_direction = new double[number_of_unknowns]();
    preconditioned_residual = new double[number_of_unknowns]();
    inverse_diagonal = new double[number_of_unknowns]();
    double right_hand_side_norm = 0.;
    out_of_iterations = false;

    for(row_index = 0; row_index < number_of_unknowns; row_index++)
        right_hand_side_norm += right_hand_side[row_index] * right_hand_side[row_index];
    right_hand_side_norm = std::sqrt(right_hand_side_norm);

    matrix.gaxpy(solution, residual_vector);
    for(row_index = 0; row_index < number_of_unknowns; row_index++)
    {
        inverse_diagonal[row_index] = 1.0 / matrix.diag[row_index].val;
        residual_vector[row_index] = right_hand_side[row_index] - residual_vector[row_index];
        preconditioned_residual[row_index] = inverse_diagonal[row_index] * residual_vector[row_index];
        search_direction[row_index] = preconditioned_residual[row_index];
    }

    residual_norm_old = dotProduct(residual_vector, preconditioned_residual);

    iteration_count = 0;
    do
    {
        matrix_times_direction = new double[number_of_unknowns]();
        matrix.gaxpy(search_direction, matrix_times_direction);

        double denominator = dotProduct(search_direction, matrix_times_direction);

        if (std::abs(denominator) < 1e-30)
        {
            std::cout << "JCG: denominator zero, stopping" << std::endl;
            break;
        }

        alpha = residual_norm_old / denominator;

        for(row_index = 0; row_index < number_of_unknowns; row_index++)
        {
            solution[row_index] += alpha * search_direction[row_index];
            residual_vector[row_index] -= alpha * matrix_times_direction[row_index];
            preconditioned_residual[row_index] = inverse_diagonal[row_index] * residual_vector[row_index];
        }

        residual_norm_new = dotProduct(residual_vector, preconditioned_residual);
        double current_residual = std::sqrt(std::abs(residual_norm_new));

        if (current_residual < solver_tolerance * right_hand_side_norm)
            break;

        beta = residual_norm_new / residual_norm_old;

        for(row_index = 0; row_index < number_of_unknowns; row_index++)
            search_direction[row_index] = preconditioned_residual[row_index] + beta * search_direction[row_index];

        residual_norm_old = residual_norm_new;

        iteration_count++;
        if (iteration_count > maximum_iterations)
        {
            out_of_iterations = true;
            break;
        }
        delete[] matrix_times_direction;
    } while (true);

    delete[] residual_vector;
    delete[] search_direction;
    delete[] matrix_times_direction;
    delete[] preconditioned_residual;
    delete[] inverse_diagonal;

    iterations = iteration_count;
    residual = std::sqrt(residual_norm_new);
}

void RealLinearSolver::solveWithSSOR(double omega)
{
    double *residual_vector = new double[number_of_unknowns]();
    double *preconditioned_residual = new double[number_of_unknowns]();
    double *search_direction = new double[number_of_unknowns]();
    double *matrix_times_direction = new double[number_of_unknowns]();

    double alpha, beta;
    double residual_norm_old, residual_norm_new;
    double right_hand_side_norm = 0.0, current_residual = 0.0;

    out_of_iterations = false;

    for (unsigned long index = 0; index < number_of_unknowns; index++)
        right_hand_side_norm += right_hand_side[index] * right_hand_side[index];
    right_hand_side_norm = std::sqrt(right_hand_side_norm);

    matrix.gaxpy(solution, residual_vector);
    for (unsigned long index = 0; index < number_of_unknowns; index++)
        residual_vector[index] = right_hand_side[index] - residual_vector[index];

    applySSORPreconditioner(residual_vector, preconditioned_residual, omega);

    for (unsigned long index = 0; index < number_of_unknowns; index++)
        search_direction[index] = preconditioned_residual[index];

    residual_norm_old = dotProduct(residual_vector, preconditioned_residual);

    iterations = 0;
    do
    {
        matrix.gaxpy(search_direction, matrix_times_direction);

        double denominator = dotProduct(search_direction, matrix_times_direction);

        if (std::abs(denominator) < 1e-30)
        {
            std::cout << "SSOR-CG: denominator zero, stopping" << std::endl;
            break;
        }

        alpha = residual_norm_old / denominator;

        for (unsigned long index = 0; index < number_of_unknowns; index++)
        {
            solution[index] += alpha * search_direction[index];
            residual_vector[index] -= alpha * matrix_times_direction[index];
        }

        applySSORPreconditioner(residual_vector, preconditioned_residual, omega);

        residual_norm_new = dotProduct(residual_vector, preconditioned_residual);
        current_residual = std::sqrt(std::abs(residual_norm_new));

        if (current_residual < solver_tolerance * right_hand_side_norm)
            break;

        beta = residual_norm_new / residual_norm_old;

        for (unsigned long index = 0; index < number_of_unknowns; index++)
            search_direction[index] = preconditioned_residual[index] + beta * search_direction[index];

        residual_norm_old = residual_norm_new;

        iterations++;
        if (iterations > maximum_iterations)
        {
            out_of_iterations = true;
            break;
        }
    } while (true);

    delete[] residual_vector;
    delete[] preconditioned_residual;
    delete[] search_direction;
    delete[] matrix_times_direction;

    residual = current_residual;
}
