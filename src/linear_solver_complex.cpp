#include <cmath>
#include <iostream>
#include "linear_solver_complex.h"

complex<double> ComplexLinearSolver::dotProduct(complex<double> *vector_one, complex<double> *vector_two)
{
    complex<double> result = 0;
    for (unsigned long index = 0; index < number_of_unknowns; index++)
        result += vector_one[index] * vector_two[index];
    return result;
}

void ComplexLinearSolver::applySSORPreconditioner(complex<double> *input_vector, complex<double> *output_vector, double omega)
{
    complex<double> scale = 1.0 / (omega * (2.0 - omega));

    // Forward sweep: solve (D + omega*L) * y = r
    // L is strictly lower triangular. For j < i, A[i][j] is stored in row j at column i
    for (unsigned long index = 0; index < number_of_unknowns; index++)
    {
        complex<double> sum = input_vector[index];
        for (unsigned long j = 0; j < index; j++)
        {
            // Search for column 'index' in row j's linked list
            ComplexEntry *element = matrix.diagonal[j].first;
            while (element && element->col < index)
                element = element->next;
            if (element && element->col == index)
                sum -= omega * element->val * output_vector[j];
        }
        output_vector[index] = sum / matrix.diagonal[index].val;
    }

    // Backward sweep: solve (D + omega*L^T) * z = y
    // L^T is upper triangular. For j > i, A[j][i] = A[i][j] stored in row i
    for (unsigned long index = number_of_unknowns; index > 0; index--)
    {
        complex<double> sum = output_vector[index-1];
        ComplexEntry *element = matrix.diagonal[index-1].first;
        while (element)
        {
            sum -= omega * element->val * output_vector[element->col];
            element = element->next;
        }
        output_vector[index-1] = sum / matrix.diagonal[index-1].val;
    }

    // Final scaling: z = z / (omega * (2 - omega))
    for (unsigned long index = 0; index < number_of_unknowns; index++)
        output_vector[index] *= scale;
}

ComplexLinearSolver::ComplexLinearSolver(unsigned long n, unsigned int max_iterations, double tolerance) : matrix(n)
{
    number_of_unknowns = n;
    solution = new complex<double>[number_of_unknowns]();
    right_hand_side = new complex<double>[number_of_unknowns]();

    solver_tolerance = tolerance;
    maximum_iterations = max_iterations;
    iterations = 0;
    residual = 0.;
    out_of_iterations = false;
}

ComplexLinearSolver::~ComplexLinearSolver()
{
    delete[] solution;
    delete[] right_hand_side;
}

void ComplexLinearSolver::setTolerance(double tolerance)
{
    solver_tolerance = tolerance;
}

void ComplexLinearSolver::setMaxIterations(unsigned int max_iterations)
{
    maximum_iterations = max_iterations;
}

unsigned int ComplexLinearSolver::getIterations()
{
    return iterations;
}

void ComplexLinearSolver::setMatrixValue(unsigned long row_index, unsigned long column_index,
                                       complex<double> value)
{
    matrix.setValue(row_index, column_index, value);
}

void ComplexLinearSolver::addToMatrixValue(unsigned long row_index, unsigned long column_index,
                                         complex<double> value)
{
    matrix.addToValue(row_index, column_index, value);
}

void ComplexLinearSolver::printResult()
{
    std::cout << " - ComplexLinearSolver Data : -" << std::endl << std::endl;
    std::cout << "[Matrix]";
    matrix.printOut();
    std::cout << "[Right Hand Side]        [Solution]" << std::endl << std::endl;
    for(unsigned long index = 0; index < number_of_unknowns; index++)
        std::cout << "[" << right_hand_side[index] << "]        [" << solution[index] << "]" << std::endl;
    std::cout << std::endl << std::endl;
}

void ComplexLinearSolver::solveWithCOCG()
{
    complex<double> *residual_vector = new complex<double>[number_of_unknowns]();
    complex<double> *search_direction = new complex<double>[number_of_unknowns]();
    complex<double> *matrix_times_direction = new complex<double>[number_of_unknowns]();

    complex<double> alpha, beta;
    complex<double> residual_norm_squared_old, residual_norm_squared_new;
    double right_hand_side_norm = 0.0, current_residual = 0.0;

    out_of_iterations = false;

    for (unsigned long index = 0; index < number_of_unknowns; index++)
        right_hand_side_norm += std::norm(right_hand_side[index]);
    right_hand_side_norm = std::sqrt(right_hand_side_norm);

    matrix.gaxpy(solution, residual_vector);
    for (unsigned long index = 0; index < number_of_unknowns; index++)
        residual_vector[index] = right_hand_side[index] - residual_vector[index];

    for (unsigned long index = 0; index < number_of_unknowns; index++)
        search_direction[index] = residual_vector[index];

    residual_norm_squared_old = dotProduct(residual_vector, residual_vector);

    iterations = 0;
    do
    {
        matrix.gaxpy(search_direction, matrix_times_direction);

        complex<double> denominator = dotProduct(search_direction, matrix_times_direction);

        if (std::abs(denominator) < 1e-30)
        {
            std::cout << "COCG: denominator zero, stopping" << std::endl;
            break;
        }

        alpha = residual_norm_squared_old / denominator;

        for (unsigned long index = 0; index < number_of_unknowns; index++)
        {
            solution[index] += alpha * search_direction[index];
            residual_vector[index] -= alpha * matrix_times_direction[index];
        }

        residual_norm_squared_new = dotProduct(residual_vector, residual_vector);
        current_residual = std::sqrt(std::abs(residual_norm_squared_new));

        if (current_residual < solver_tolerance * right_hand_side_norm)
            break;

        beta = residual_norm_squared_new / residual_norm_squared_old;

        for (unsigned long index = 0; index < number_of_unknowns; index++)
            search_direction[index] = residual_vector[index] + beta * search_direction[index];

        residual_norm_squared_old = residual_norm_squared_new;

        iterations++;
        if (iterations > maximum_iterations)
        {
            out_of_iterations = true;
            break;
        }
    } while (true);

    delete[] residual_vector;
    delete[] search_direction;
    delete[] matrix_times_direction;

    residual = current_residual;
}

void ComplexLinearSolver::solveWithCOCGAndSSOR(double omega)
{
    complex<double> *residual_vector = new complex<double>[number_of_unknowns]();
    complex<double> *preconditioned_residual = new complex<double>[number_of_unknowns]();
    complex<double> *search_direction = new complex<double>[number_of_unknowns]();
    complex<double> *matrix_times_direction = new complex<double>[number_of_unknowns]();

    complex<double> alpha, beta;
    complex<double> residual_norm_squared_old, residual_norm_squared_new;
    double right_hand_side_norm = 0.0, current_residual = 0.0;

    out_of_iterations = false;

    for (unsigned long index = 0; index < number_of_unknowns; index++)
        right_hand_side_norm += std::norm(right_hand_side[index]);
    right_hand_side_norm = std::sqrt(right_hand_side_norm);

    matrix.gaxpy(solution, residual_vector);
    for (unsigned long index = 0; index < number_of_unknowns; index++)
        residual_vector[index] = right_hand_side[index] - residual_vector[index];

    applySSORPreconditioner(residual_vector, preconditioned_residual, omega);

    for (unsigned long index = 0; index < number_of_unknowns; index++)
        search_direction[index] = preconditioned_residual[index];

    residual_norm_squared_old = dotProduct(residual_vector, preconditioned_residual);

    iterations = 0;
    do
    {
        matrix.gaxpy(search_direction, matrix_times_direction);

        complex<double> denominator = dotProduct(search_direction, matrix_times_direction);

        if (std::abs(denominator) < 1e-30)
        {
            std::cout << "COCG-SSOR: denominator zero, stopping" << std::endl;
            break;
        }

        alpha = residual_norm_squared_old / denominator;

        for (unsigned long index = 0; index < number_of_unknowns; index++)
        {
            solution[index] += alpha * search_direction[index];
            residual_vector[index] -= alpha * matrix_times_direction[index];
        }

        applySSORPreconditioner(residual_vector, preconditioned_residual, omega);

        residual_norm_squared_new = dotProduct(residual_vector, preconditioned_residual);
        current_residual = std::sqrt(std::abs(residual_norm_squared_new));

        if (current_residual < solver_tolerance * right_hand_side_norm)
            break;

        beta = residual_norm_squared_new / residual_norm_squared_old;

        for (unsigned long index = 0; index < number_of_unknowns; index++)
            search_direction[index] = preconditioned_residual[index] + beta * search_direction[index];

        residual_norm_squared_old = residual_norm_squared_new;

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
