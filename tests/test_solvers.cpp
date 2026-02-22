#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <cmath>
#include <vector>
#include <complex>

#include <doctest.h>

#include "linear_solver_real.h"
#include "linear_solver_complex.h"

TEST_CASE("Sparse Matrix - Basic Operations")
{
    SUBCASE("Set and get diagonal values")
    {
        Sparse A(3);
        A.setValue(0, 0, 1.0, true);
        A.setValue(1, 1, 2.0, true);
        A.setValue(2, 2, 3.0, true);

        REQUIRE(A.diag[0].val == 1.0);
        REQUIRE(A.diag[1].val == 2.0);
        REQUIRE(A.diag[2].val == 3.0);
    }

    SUBCASE("Set off-diagonal values (symmetric)")
    {
        Sparse A(3);
        A.setValue(0, 0, 1.0, true);
        A.setValue(1, 1, 2.0, true);
        A.setValue(2, 2, 1.0, true);
        A.setValue(0, 1, 5.0, true);
        A.setValue(0, 2, -3.0, true);
        A.setValue(1, 2, 2.0, true);

        double x[3] = {1.0, 1.0, 1.0};
        double y[3] = {0.0};

        A.gaxpy(x, y);

        // y[0] = 1*1 + 5*1 + (-3)*1 = 1 + 5 - 3 = 3
        // y[1] = 5*1 + 2*1 + 2*1 = 5 + 2 + 2 = 9
        // y[2] = (-3)*1 + 2*1 + 1*1 = -3 + 2 + 1 = 0
        REQUIRE(y[0] == doctest::Approx(3.0));
        REQUIRE(y[1] == doctest::Approx(9.0));
        REQUIRE(y[2] == doctest::Approx(0.0));
    }

    SUBCASE("gaxpy computes A*x correctly")
    {
        Sparse A(2);
        A.setValue(0, 0, 4.0, true);
        A.setValue(1, 1, 9.0, true);
        A.setValue(0, 1, 1.0, true);

        double x[2] = {2.0, 3.0};
        double y[2] = {0.0};

        A.gaxpy(x, y);

        double expected0 = 4.0 * 2.0 + 1.0 * 3.0; // 8 + 3 = 11
        double expected1 = 1.0 * 2.0 + 9.0 * 3.0; // 2 + 27 = 29

        REQUIRE(y[0] == doctest::Approx(expected0));
        REQUIRE(y[1] == doctest::Approx(expected1));
    }
}

TEST_CASE("CG Solver - Conjugate Gradient")
{
    SUBCASE("3x3 system solves correctly")
    {
        unsigned long n = 3;
        RealLinearSolver prob(n);
        prob.setMaxIterations(10000);
        prob.setTolerance(1e-15);

        // Matrix: symmetric positive definite
        // |  1  5 -3 |
        // |  5  2  2 |
        // | -3  2  1 |
        prob.setMatrixValue(0, 0, 1.0, true);
        prob.setMatrixValue(1, 1, 2.0, true);
        prob.setMatrixValue(2, 2, 1.0, true);
        prob.setMatrixValue(0, 1, 5.0, true);
        prob.setMatrixValue(0, 2, -3.0, true);
        prob.setMatrixValue(1, 2, 2.0, true);

        // RHS: B = [8, -2, 4]
        prob.right_hand_side[0] = 8.0;
        prob.right_hand_side[1] = -2.0;
        prob.right_hand_side[2] = 4.0;

        prob.solveWithCG();

        // Verify solution: Ax = b
        double residual = 0.0;
        double b_norm = 0.0;
        for (unsigned long index = 0; index < n; index++)
        {
            b_norm += prob.right_hand_side[index] * prob.right_hand_side[index];
        }
        b_norm = std::sqrt(b_norm);

        // Compute Ax - b
        double Ax[3] = {0.0};
        prob.matrix.gaxpy(prob.solution, Ax);
        for (unsigned long index = 0; index < n; index++)
        {
            double ri = Ax[index] - prob.right_hand_side[index];
            residual += ri * ri;
        }
        residual = std::sqrt(residual);

        double relative_error = residual / b_norm;

        DOCTEST_INFO("Solution: X = [" << prob.solution[0] << ", " << prob.solution[1] << ", " << prob.solution[2] << "]");
        DOCTEST_INFO("Relative residual: " << relative_error);

        REQUIRE(relative_error < 1e-10);
    }

    SUBCASE("JCG solver converges")
    {
        unsigned long n = 3;
        RealLinearSolver prob(n);
        prob.setMaxIterations(10000);
        prob.setTolerance(1e-15);

        prob.setMatrixValue(0, 0, 1.0, true);
        prob.setMatrixValue(1, 1, 2.0, true);
        prob.setMatrixValue(2, 2, 1.0, true);
        prob.setMatrixValue(0, 1, 5.0, true);
        prob.setMatrixValue(0, 2, -3.0, true);
        prob.setMatrixValue(1, 2, 2.0, true);

        prob.right_hand_side[0] = 8.0;
        prob.right_hand_side[1] = -2.0;
        prob.right_hand_side[2] = 4.0;

        prob.solveWithJCG();

        double residual = 0.0;
        double Ax[3] = {0.0};
        prob.matrix.gaxpy(prob.solution, Ax);
        for (unsigned long index = 0; index < n; index++)
        {
            double ri = Ax[index] - prob.right_hand_side[index];
            residual += ri * ri;
        }
        residual = std::sqrt(residual);

        REQUIRE(residual < 1e-8);

        REQUIRE(prob.solution[0] == doctest::Approx(-2.0 / 3.0));
        REQUIRE(prob.solution[1] == doctest::Approx(4.0 / 3.0));
        REQUIRE(prob.solution[2] == doctest::Approx(-2.0 / 3.0));
    }

    SUBCASE("SSOR-CG solver converges with omega=1.0")
    {
        unsigned long n = 3;
        RealLinearSolver prob(n);
        prob.setMaxIterations(10000);
        prob.setTolerance(1e-10);

        prob.setMatrixValue(0, 0, 1.0, true);
        prob.setMatrixValue(1, 1, 2.0, true);
        prob.setMatrixValue(2, 2, 1.0, true);
        prob.setMatrixValue(0, 1, 5.0, true);
        prob.setMatrixValue(0, 2, -3.0, true);
        prob.setMatrixValue(1, 2, 2.0, true);

        prob.right_hand_side[0] = 8.0;
        prob.right_hand_side[1] = -2.0;
        prob.right_hand_side[2] = 4.0;

        prob.solveWithSSOR(1.0);

        double residual = 0.0;
        double Ax[3] = {0.0};
        prob.matrix.gaxpy(prob.solution, Ax);
        for (unsigned long index = 0; index < n; index++)
        {
            double ri = Ax[index] - prob.right_hand_side[index];
            residual += ri * ri;
        }
        residual = std::sqrt(residual);

        DOCTEST_INFO("Solution: X = [" << prob.solution[0] << ", " << prob.solution[1] << ", " << prob.solution[2] << "]");
        DOCTEST_INFO("Residual: " << residual);

        REQUIRE(residual < 1e-6);
    }

    SUBCASE("SSOR-CG solver converges with omega=1.2 (over-relaxation)")
    {
        unsigned long n = 3;
        RealLinearSolver prob(n);
        prob.setMaxIterations(10000);
        prob.setTolerance(1e-10);

        prob.setMatrixValue(0, 0, 1.0, true);
        prob.setMatrixValue(1, 1, 2.0, true);
        prob.setMatrixValue(2, 2, 1.0, true);
        prob.setMatrixValue(0, 1, 5.0, true);
        prob.setMatrixValue(0, 2, -3.0, true);
        prob.setMatrixValue(1, 2, 2.0, true);

        prob.right_hand_side[0] = 8.0;
        prob.right_hand_side[1] = -2.0;
        prob.right_hand_side[2] = 4.0;

        prob.solveWithSSOR(1.2);

        double residual = 0.0;
        double Ax[3] = {0.0};
        prob.matrix.gaxpy(prob.solution, Ax);
        for (unsigned long index = 0; index < n; index++)
        {
            double ri = Ax[index] - prob.right_hand_side[index];
            residual += ri * ri;
        }
        residual = std::sqrt(residual);

        DOCTEST_INFO("Solution: X = [" << prob.solution[0] << ", " << prob.solution[1] << ", " << prob.solution[2] << "]");
        DOCTEST_INFO("Residual: " << residual);
        DOCTEST_INFO("Iterations: " << prob.getIterations());

        REQUIRE(residual < 1e-6);
        REQUIRE(prob.solution[0] == doctest::Approx(-2.0 / 3.0).epsilon(1e-6));
        REQUIRE(prob.solution[1] == doctest::Approx(4.0 / 3.0).epsilon(1e-6));
        REQUIRE(prob.solution[2] == doctest::Approx(-2.0 / 3.0).epsilon(1e-6));
    }

    SUBCASE("SSOR-CG solver converges with omega=0.8 (under-relaxation)")
    {
        unsigned long n = 3;
        RealLinearSolver prob(n);
        prob.setMaxIterations(10000);
        prob.setTolerance(1e-10);

        prob.setMatrixValue(0, 0, 1.0, true);
        prob.setMatrixValue(1, 1, 2.0, true);
        prob.setMatrixValue(2, 2, 1.0, true);
        prob.setMatrixValue(0, 1, 5.0, true);
        prob.setMatrixValue(0, 2, -3.0, true);
        prob.setMatrixValue(1, 2, 2.0, true);

        prob.right_hand_side[0] = 8.0;
        prob.right_hand_side[1] = -2.0;
        prob.right_hand_side[2] = 4.0;

        prob.solveWithSSOR(0.8);

        double residual = 0.0;
        double Ax[3] = {0.0};
        prob.matrix.gaxpy(prob.solution, Ax);
        for (unsigned long index = 0; index < n; index++)
        {
            double ri = Ax[index] - prob.right_hand_side[index];
            residual += ri * ri;
        }
        residual = std::sqrt(residual);

        DOCTEST_INFO("Solution: X = [" << prob.solution[0] << ", " << prob.solution[1] << ", " << prob.solution[2] << "]");
        DOCTEST_INFO("Residual: " << residual);
        DOCTEST_INFO("Iterations: " << prob.getIterations());

        REQUIRE(residual < 1e-6);
    }
}

TEST_CASE("Matrix - Replace vs Add")
{
    SUBCASE("replace=true overwrites value")
    {
        Sparse A(2);
        A.setValue(0, 0, 5.0, true); // replace
        A.setValue(0, 0, 3.0, true); // replace again

        REQUIRE(A.diag[0].val == 3.0);
    }

    SUBCASE("replace=false adds to value")
    {
        Sparse A(2);
        A.setValue(0, 0, 5.0, true);
        A.setValue(0, 0, 3.0, false); // add

        REQUIRE(A.diag[0].val == 8.0);
    }
}

TEST_CASE("Complex Sparse Matrix - Basic Operations")
{
    SUBCASE("Set and get diagonal values")
    {
        ComplexSparseMatrix A(3);
        A.setValue(0, 0, complex<double>(1.0, 0.0), true);
        A.setValue(1, 1, complex<double>(2.0, 0.0), true);
        A.setValue(2, 2, complex<double>(3.0, 0.0), true);

        REQUIRE(A.diagonal[0].val.real() == doctest::Approx(1.0));
        REQUIRE(A.diagonal[1].val.real() == doctest::Approx(2.0));
        REQUIRE(A.diagonal[2].val.real() == doctest::Approx(3.0));
    }

    SUBCASE("Set symmetric off-diagonal values (for eddy current)")
    {
        ComplexSparseMatrix A(2);
        A.setValue(0, 0, complex<double>(1.0, 0.0), true);
        A.setValue(1, 1, complex<double>(2.0, 0.0), true);
        A.setValue(0, 1, complex<double>(3.0, 4.0), true);

        complex<double> x[2] = {complex<double>(1.0, 0.0), complex<double>(1.0, 0.0)};
        complex<double> y[2] = {0.0, 0.0};

        A.gaxpy(x, y);

        // y[0] = 1*1 + (3+4i)*1 = 1 + 3 + 4i = 4 + 4i
        // y[1] = (3+4i)*1 + 2*1 = 3+4i + 2 = 5 + 4i
        REQUIRE(y[0].real() == doctest::Approx(4.0));
        REQUIRE(y[0].imag() == doctest::Approx(4.0));
        REQUIRE(y[1].real() == doctest::Approx(5.0));
        REQUIRE(y[1].imag() == doctest::Approx(4.0));
    }
}

TEST_CASE("COCG Solver - Conjugate Orthogonal CG for Symmetric Complex")
{
    SUBCASE("COCG solves 2x2 symmetric complex system (eddy current)")
    {
        unsigned long n = 2;
        ComplexLinearSolver prob(n);
        prob.setMaxIterations(1000);
        prob.setTolerance(1e-10);

        // Symmetric complex matrix (typical for eddy current problems)
        // | 3+0i   1+2i |
        // | 1+2i   4+0i |
        prob.setMatrixValue(0, 0, complex<double>(3.0, 0.0), true);
        prob.setMatrixValue(1, 1, complex<double>(4.0, 0.0), true);
        prob.setMatrixValue(0, 1, complex<double>(1.0, 2.0), true);

        prob.right_hand_side[0] = complex<double>(10.0, 0.0);
        prob.right_hand_side[1] = complex<double>(8.0, 0.0);

        prob.solveWithCOCG();

        complex<double> Ax[2] = {0.0, 0.0};
        prob.matrix.gaxpy(prob.solution, Ax);

        double residual = 0.0;
        for (unsigned long i = 0; i < n; i++)
        {
            complex<double> r = Ax[i] - prob.right_hand_side[i];
            residual += std::norm(r);
        }
        residual = std::sqrt(residual);

        DOCTEST_INFO("Solution: X = [" << prob.solution[0] << ", " << prob.solution[1] << "]");
        DOCTEST_INFO("Residual: " << residual);
        DOCTEST_INFO("Iterations: " << prob.getIterations());

        REQUIRE(residual < 1e-6);
    }

    SUBCASE("COCG solves 3x3 symmetric complex system")
    {
        unsigned long n = 3;
        ComplexLinearSolver prob(n);
        prob.setMaxIterations(1000);
        prob.setTolerance(1e-12);

        // Symmetric complex matrix
        // | 4+0i   1+1i   0+0i |
        // | 1+1i   3+0i   2-1i |
        // | 0+0i   2-1i   5+0i |
        prob.setMatrixValue(0, 0, complex<double>(4.0, 0.0), true);
        prob.setMatrixValue(1, 1, complex<double>(3.0, 0.0), true);
        prob.setMatrixValue(2, 2, complex<double>(5.0, 0.0), true);
        prob.setMatrixValue(0, 1, complex<double>(1.0, 1.0), true);
        prob.setMatrixValue(1, 2, complex<double>(2.0, -1.0), true);

        prob.right_hand_side[0] = complex<double>(5.0, 1.0);
        prob.right_hand_side[1] = complex<double>(6.0, 2.0);
        prob.right_hand_side[2] = complex<double>(7.0, -1.0);

        prob.solveWithCOCG();

        complex<double> Ax[3] = {0.0, 0.0, 0.0};
        prob.matrix.gaxpy(prob.solution, Ax);

        double residual = 0.0;
        for (unsigned long i = 0; i < n; i++)
        {
            complex<double> r = Ax[i] - prob.right_hand_side[i];
            residual += std::norm(r);
        }
        residual = std::sqrt(residual);

        DOCTEST_INFO("Solution: X = [" << prob.solution[0] << ", " << prob.solution[1] << ", " << prob.solution[2] << "]");
        DOCTEST_INFO("Residual: " << residual);

        REQUIRE(residual < 1e-8);
    }

    SUBCASE("COCG-SSOR solves 3x3 symmetric complex system faster")
    {
        unsigned long n = 3;
        ComplexLinearSolver prob(n);
        prob.setMaxIterations(1000);
        prob.setTolerance(1e-12);

        prob.setMatrixValue(0, 0, complex<double>(4.0, 0.0), true);
        prob.setMatrixValue(1, 1, complex<double>(3.0, 0.0), true);
        prob.setMatrixValue(2, 2, complex<double>(5.0, 0.0), true);
        prob.setMatrixValue(0, 1, complex<double>(1.0, 1.0), true);
        prob.setMatrixValue(1, 2, complex<double>(2.0, -1.0), true);

        prob.right_hand_side[0] = complex<double>(5.0, 1.0);
        prob.right_hand_side[1] = complex<double>(6.0, 2.0);
        prob.right_hand_side[2] = complex<double>(7.0, -1.0);

        prob.solveWithCOCGAndSSOR(1.0);

        complex<double> Ax[3] = {0.0, 0.0, 0.0};
        prob.matrix.gaxpy(prob.solution, Ax);

        double residual = 0.0;
        for (unsigned long i = 0; i < n; i++)
        {
            complex<double> r = Ax[i] - prob.right_hand_side[i];
            residual += std::norm(r);
        }
        residual = std::sqrt(residual);

        DOCTEST_INFO("Solution: X = [" << prob.solution[0] << ", " << prob.solution[1] << ", " << prob.solution[2] << "]");
        DOCTEST_INFO("Residual: " << residual);
        DOCTEST_INFO("Iterations: " << prob.getIterations());

        REQUIRE(residual < 1e-8);
    }

    SUBCASE("COCG-SSOR solves 3x3 symmetric complex system with omega=1.2")
    {
        unsigned long n = 3;
        ComplexLinearSolver prob(n);
        prob.setMaxIterations(1000);
        prob.setTolerance(1e-12);

        prob.setMatrixValue(0, 0, complex<double>(4.0, 0.0), true);
        prob.setMatrixValue(1, 1, complex<double>(3.0, 0.0), true);
        prob.setMatrixValue(2, 2, complex<double>(5.0, 0.0), true);
        prob.setMatrixValue(0, 1, complex<double>(1.0, 1.0), true);
        prob.setMatrixValue(1, 2, complex<double>(2.0, -1.0), true);

        prob.right_hand_side[0] = complex<double>(5.0, 1.0);
        prob.right_hand_side[1] = complex<double>(6.0, 2.0);
        prob.right_hand_side[2] = complex<double>(7.0, -1.0);

        prob.solveWithCOCGAndSSOR(1.2);

        complex<double> Ax[3] = {0.0, 0.0, 0.0};
        prob.matrix.gaxpy(prob.solution, Ax);

        double residual = 0.0;
        for (unsigned long i = 0; i < n; i++)
        {
            complex<double> r = Ax[i] - prob.right_hand_side[i];
            residual += std::norm(r);
        }
        residual = std::sqrt(residual);

        DOCTEST_INFO("Solution: X = [" << prob.solution[0] << ", " << prob.solution[1] << ", " << prob.solution[2] << "]");
        DOCTEST_INFO("Residual: " << residual);
        DOCTEST_INFO("Iterations: " << prob.getIterations());

        REQUIRE(residual < 1e-8);
    }
}
