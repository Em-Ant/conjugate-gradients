#ifndef LINEAR_SOLVER_H_INCLUDED
#define LINEAR_SOLVER_H_INCLUDED

#include "sparse_matrix.h"
#include <iostream>
#include <cmath>

template<typename Traits>
class LinearSolver {
public:
    using Scalar = typename Traits::Scalar;
    using Matrix = SparseMatrix<Traits>;
    
    unsigned long number_of_unknowns;
    Scalar *solution;
    Scalar *right_hand_side;
    Matrix matrix;
    
    LinearSolver(unsigned long n, unsigned int max_iterations = 50000,
                 double tolerance = 1e-15, unsigned int iteration_update = 2000)
        : number_of_unknowns(n), solution(nullptr), right_hand_side(nullptr),
          matrix(n),
          maximum_iterations(max_iterations), solver_tolerance(tolerance),
          iteration_update_frequency(iteration_update)
    {
        solution = new Scalar[number_of_unknowns]();
        right_hand_side = new Scalar[number_of_unknowns]();
        iterations = 0;
        residual = 0.0;
        out_of_iterations = false;
    }
    
    ~LinearSolver() {
        delete[] solution;
        delete[] right_hand_side;
    }
    
    void setTolerance(double tolerance) {
        solver_tolerance = tolerance;
    }
    
    void setMaxIterations(unsigned int max_iterations) {
        maximum_iterations = max_iterations;
    }
    
    void setIterationUpdate(unsigned int iteration_update) {
        iteration_update_frequency = iteration_update;
    }
    
    void setMatrixValue(unsigned long row_index, unsigned long column_index, const Scalar& value) {
        matrix.setValue(row_index, column_index, value);
    }
    
    void addToMatrixValue(unsigned long row_index, unsigned long column_index, const Scalar& value) {
        matrix.addToValue(row_index, column_index, value);
    }
    
    void printResult() {
        std::cout << " - Linear Solver Data : -" << std::endl << std::endl;
        std::cout << "[Matrix]";
        matrix.printOut();
        std::cout << "[Right Hand Side]        [Solution]" << std::endl << std::endl;
        for (unsigned long index = 0; index < number_of_unknowns; index++) {
            std::cout << "[" << right_hand_side[index] << "]        [" << solution[index] << "]" << std::endl;
        }
        std::cout << std::endl << std::endl;
    }
    
    unsigned int getIterations() {
        return iterations;
    }
    
    Scalar dotProduct(Scalar *vector_one, Scalar *vector_two) {
        Scalar result = Traits::zero();
        for (unsigned long index = 0; index < number_of_unknowns; index++) {
            result += Traits::conjugate(vector_one[index]) * vector_two[index];
        }
        return result;
    }
    
    void solveWithCG() {
        static_assert(Traits::supports_symmetric || Traits::supports_hermitian,
                      "CG requires symmetric or Hermitian positive definite matrix");
        
        unsigned long row_index, iteration_count;
        Scalar *matrix_times_direction;
        Scalar *search_direction;
        Scalar *residual_vector;
        Scalar alpha, beta;
        Scalar residual_norm_old, residual_norm_new = Traits::zero();
        
        residual_vector = new Scalar[number_of_unknowns]();
        search_direction = new Scalar[number_of_unknowns]();
        matrix_times_direction = new Scalar[number_of_unknowns]();
        
        double right_hand_side_norm = 0.0;
        out_of_iterations = false;
        
        for (row_index = 0; row_index < number_of_unknowns; row_index++) {
            right_hand_side_norm += Traits::norm(right_hand_side[row_index]);
        }
        right_hand_side_norm = std::sqrt(right_hand_side_norm);
        
        matrix.gaxpy(solution, residual_vector);
        for (row_index = 0; row_index < number_of_unknowns; row_index++) {
            residual_vector[row_index] = right_hand_side[row_index] - residual_vector[row_index];
            search_direction[row_index] = residual_vector[row_index];
        }
        
        residual_norm_old = dotProduct(residual_vector, residual_vector);
        
        iteration_count = 0;
        do {
            for (unsigned long i = 0; i < number_of_unknowns; i++) {
                matrix_times_direction[i] = Traits::zero();
            }
            matrix.gaxpy(search_direction, matrix_times_direction);
            
            Scalar denominator = dotProduct(search_direction, matrix_times_direction);
            
            if (Traits::abs(denominator) < 1e-30) {
                std::cout << "CG: denominator zero, stopping" << std::endl;
                break;
            }
            
            alpha = residual_norm_old / denominator;
            
            for (row_index = 0; row_index < number_of_unknowns; row_index++) {
                solution[row_index] += alpha * search_direction[row_index];
                residual_vector[row_index] -= alpha * matrix_times_direction[row_index];
            }
            
            residual_norm_new = dotProduct(residual_vector, residual_vector);
            double current_residual = std::sqrt(Traits::abs(residual_norm_new));
            
            if (current_residual < solver_tolerance * right_hand_side_norm) {
                break;
            }
            
            beta = residual_norm_new / residual_norm_old;
            
            for (row_index = 0; row_index < number_of_unknowns; row_index++) {
                search_direction[row_index] = residual_vector[row_index] + beta * search_direction[row_index];
            }
            
            residual_norm_old = residual_norm_new;
            
            iteration_count++;
            if (iteration_count > maximum_iterations) {
                out_of_iterations = true;
                break;
            }
        } while (true);
        
        delete[] residual_vector;
        delete[] search_direction;
        delete[] matrix_times_direction;
        
        iterations = iteration_count;
        residual = std::sqrt(Traits::abs(residual_norm_new));
    }
    
    void solveWithCOCG() {
        static_assert(Traits::is_complex && Traits::supports_symmetric,
                      "COCG requires complex symmetric matrix");
        
        Scalar *residual_vector = new Scalar[number_of_unknowns]();
        Scalar *search_direction = new Scalar[number_of_unknowns]();
        Scalar *matrix_times_direction = new Scalar[number_of_unknowns]();
        
        Scalar alpha, beta;
        Scalar residual_norm_squared_old, residual_norm_squared_new;
        double right_hand_side_norm = 0.0, current_residual = 0.0;
        
        out_of_iterations = false;
        
        for (unsigned long index = 0; index < number_of_unknowns; index++) {
            right_hand_side_norm += Traits::norm(right_hand_side[index]);
        }
        right_hand_side_norm = std::sqrt(right_hand_side_norm);
        
        matrix.gaxpy(solution, residual_vector);
        for (unsigned long index = 0; index < number_of_unknowns; index++) {
            residual_vector[index] = right_hand_side[index] - residual_vector[index];
        }
        
        for (unsigned long index = 0; index < number_of_unknowns; index++) {
            search_direction[index] = residual_vector[index];
        }
        
        residual_norm_squared_old = dotProduct(residual_vector, residual_vector);
        
        iterations = 0;
        do {
            for (unsigned long i = 0; i < number_of_unknowns; i++) {
                matrix_times_direction[i] = Traits::zero();
            }
            matrix.gaxpy(search_direction, matrix_times_direction);
            
            Scalar denominator = dotProduct(search_direction, matrix_times_direction);
            
            if (Traits::abs(denominator) < 1e-30) {
                std::cout << "COCG: denominator zero, stopping" << std::endl;
                break;
            }
            
            alpha = residual_norm_squared_old / denominator;
            
            for (unsigned long index = 0; index < number_of_unknowns; index++) {
                solution[index] += alpha * search_direction[index];
                residual_vector[index] -= alpha * matrix_times_direction[index];
            }
            
            residual_norm_squared_new = dotProduct(residual_vector, residual_vector);
            current_residual = std::sqrt(Traits::abs(residual_norm_squared_new));
            
            if (current_residual < solver_tolerance * right_hand_side_norm) {
                break;
            }
            
            beta = residual_norm_squared_new / residual_norm_squared_old;
            
            for (unsigned long index = 0; index < number_of_unknowns; index++) {
                search_direction[index] = residual_vector[index] + beta * search_direction[index];
            }
            
            residual_norm_squared_old = residual_norm_squared_new;
            
            iterations++;
            if (iterations > maximum_iterations) {
                out_of_iterations = true;
                break;
            }
        } while (true);
        
        delete[] residual_vector;
        delete[] search_direction;
        delete[] matrix_times_direction;
        
        residual = current_residual;
    }
    
    void solveWithJCG() {
        static_assert(!Traits::is_complex,
                      "JCG requires real symmetric matrix");
        
        unsigned long row_index, iteration_count;
        Scalar *matrix_times_direction;
        Scalar *search_direction;
        Scalar *residual_vector;
        Scalar *preconditioned_residual;
        Scalar alpha, beta;
        Scalar residual_norm_old, residual_norm_new = Traits::zero();
        Scalar *inverse_diagonal;
        
        residual_vector = new Scalar[number_of_unknowns]();
        search_direction = new Scalar[number_of_unknowns]();
        matrix_times_direction = new Scalar[number_of_unknowns]();
        preconditioned_residual = new Scalar[number_of_unknowns]();
        inverse_diagonal = new Scalar[number_of_unknowns]();
        
        double right_hand_side_norm = 0.0;
        out_of_iterations = false;
        
        for (row_index = 0; row_index < number_of_unknowns; row_index++) {
            right_hand_side_norm += Traits::norm(right_hand_side[row_index]);
        }
        right_hand_side_norm = std::sqrt(right_hand_side_norm);
        
        matrix.gaxpy(solution, residual_vector);
        for (row_index = 0; row_index < number_of_unknowns; row_index++) {
            inverse_diagonal[row_index] = 1.0 / matrix.diag[row_index].val;
            residual_vector[row_index] = right_hand_side[row_index] - residual_vector[row_index];
            preconditioned_residual[row_index] = inverse_diagonal[row_index] * residual_vector[row_index];
            search_direction[row_index] = preconditioned_residual[row_index];
        }
        
        residual_norm_old = dotProduct(residual_vector, preconditioned_residual);
        
        iteration_count = 0;
        do {
            for (unsigned long i = 0; i < number_of_unknowns; i++) {
                matrix_times_direction[i] = Traits::zero();
            }
            matrix.gaxpy(search_direction, matrix_times_direction);
            
            Scalar denominator = dotProduct(search_direction, matrix_times_direction);
            
            if (Traits::abs(denominator) < 1e-30) {
                std::cout << "JCG: denominator zero, stopping" << std::endl;
                break;
            }
            
            alpha = residual_norm_old / denominator;
            
            for (row_index = 0; row_index < number_of_unknowns; row_index++) {
                solution[row_index] += alpha * search_direction[row_index];
                residual_vector[row_index] -= alpha * matrix_times_direction[row_index];
                preconditioned_residual[row_index] = inverse_diagonal[row_index] * residual_vector[row_index];
            }
            
            residual_norm_new = dotProduct(residual_vector, preconditioned_residual);
            double current_residual = std::sqrt(Traits::abs(residual_norm_new));
            
            if (current_residual < solver_tolerance * right_hand_side_norm) {
                break;
            }
            
            beta = residual_norm_new / residual_norm_old;
            
            for (row_index = 0; row_index < number_of_unknowns; row_index++) {
                search_direction[row_index] = preconditioned_residual[row_index] + beta * search_direction[row_index];
            }
            
            residual_norm_old = residual_norm_new;
            
            iteration_count++;
            if (iteration_count > maximum_iterations) {
                out_of_iterations = true;
                break;
            }
        } while (true);
        
        delete[] residual_vector;
        delete[] search_direction;
        delete[] matrix_times_direction;
        delete[] preconditioned_residual;
        delete[] inverse_diagonal;
        
        iterations = iteration_count;
        residual = std::sqrt(Traits::abs(residual_norm_new));
    }
    
    void solveWithSSOR(double omega) {
        Scalar *residual_vector = new Scalar[number_of_unknowns]();
        Scalar *preconditioned_residual = new Scalar[number_of_unknowns]();
        Scalar *search_direction = new Scalar[number_of_unknowns]();
        Scalar *matrix_times_direction = new Scalar[number_of_unknowns]();
        
        Scalar alpha, beta;
        Scalar residual_norm_old, residual_norm_new = Traits::zero();
        double right_hand_side_norm = 0.0, current_residual = 0.0;
        
        out_of_iterations = false;
        
        for (unsigned long index = 0; index < number_of_unknowns; index++) {
            right_hand_side_norm += Traits::norm(right_hand_side[index]);
        }
        right_hand_side_norm = std::sqrt(right_hand_side_norm);
        
        matrix.gaxpy(solution, residual_vector);
        for (unsigned long index = 0; index < number_of_unknowns; index++) {
            residual_vector[index] = right_hand_side[index] - residual_vector[index];
        }
        
        applySSORPreconditioner(residual_vector, preconditioned_residual, omega);
        
        for (unsigned long index = 0; index < number_of_unknowns; index++) {
            search_direction[index] = preconditioned_residual[index];
        }
        
        residual_norm_old = dotProduct(residual_vector, preconditioned_residual);
        
        iterations = 0;
        do {
            for (unsigned long i = 0; i < number_of_unknowns; i++) {
                matrix_times_direction[i] = Traits::zero();
            }
            matrix.gaxpy(search_direction, matrix_times_direction);
            
            Scalar denominator = dotProduct(search_direction, matrix_times_direction);
            
            if (Traits::abs(denominator) < 1e-30) {
                std::cout << "SSOR-CG: denominator zero, stopping" << std::endl;
                break;
            }
            
            alpha = residual_norm_old / denominator;
            
            for (unsigned long index = 0; index < number_of_unknowns; index++) {
                solution[index] += alpha * search_direction[index];
                residual_vector[index] -= alpha * matrix_times_direction[index];
            }
            
            applySSORPreconditioner(residual_vector, preconditioned_residual, omega);
            
            residual_norm_new = dotProduct(residual_vector, preconditioned_residual);
            current_residual = std::sqrt(Traits::abs(residual_norm_new));
            
            if (current_residual < solver_tolerance * right_hand_side_norm) {
                break;
            }
            
            beta = residual_norm_new / residual_norm_old;
            
            for (unsigned long index = 0; index < number_of_unknowns; index++) {
                search_direction[index] = preconditioned_residual[index] + beta * search_direction[index];
            }
            
            residual_norm_old = residual_norm_new;
            
            iterations++;
            if (iterations > maximum_iterations) {
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

private:
    void applySSORPreconditioner(Scalar *input_vector, Scalar *output_vector, double omega) {
        unsigned long index;
        Scalar scale = 1.0 / (omega * (2.0 - omega));
        typename Matrix::Entry *element;
        
        for (index = 0; index < number_of_unknowns; index++) {
            output_vector[index] = Traits::zero();
        }
        
        for (index = 0; index < number_of_unknowns; index++) {
            Scalar sum = input_vector[index];
            for (unsigned long j = 0; j < index; j++) {
                element = matrix.diag[j].first;
                while (element && element->col < index) {
                    element = element->next;
                }
                if (element && element->col == index) {
                    sum -= omega * element->val * output_vector[j];
                }
            }
            output_vector[index] = sum / matrix.diag[index].val;
        }
        
        for (index = number_of_unknowns; index > 0; index--) {
            Scalar sum = output_vector[index - 1];
            element = matrix.diag[index - 1].first;
            while (element) {
                sum -= omega * element->val * output_vector[element->col];
                element = element->next;
            }
            output_vector[index - 1] = sum / matrix.diag[index - 1].val;
        }
        
        for (index = 0; index < number_of_unknowns; index++) {
            output_vector[index] *= scale;
        }
    }
    
    unsigned int maximum_iterations;
    double solver_tolerance;
    unsigned int iteration_update_frequency;
    unsigned int iterations;
    double residual;
    bool out_of_iterations;
};

using RealLinearSolver = LinearSolver<RealTraits>;
using ComplexLinearSolver = LinearSolver<ComplexSymmetricTraits>;

#endif
