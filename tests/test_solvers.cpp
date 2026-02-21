#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <cmath>
#include <vector>
#include "symsparseprob.h"

using Catch::Approx;

TEST_CASE("Sparse Matrix - Basic Operations", "[sparse]")
{
    SECTION("Set and get diagonal values")
    {
        Sparse A(3);
        A.setValue(0, 0, 1.0, true);
        A.setValue(1, 1, 2.0, true);
        A.setValue(2, 2, 3.0, true);

        REQUIRE(A.diag[0].val == 1.0);
        REQUIRE(A.diag[1].val == 2.0);
        REQUIRE(A.diag[2].val == 3.0);
    }

    SECTION("Set off-diagonal values (symmetric)")
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
        REQUIRE(y[0] == Approx(3.0));
        REQUIRE(y[1] == Approx(9.0));
        REQUIRE(y[2] == Approx(0.0));
    }

    SECTION("gaxpy computes A*x correctly")
    {
        Sparse A(2);
        A.setValue(0, 0, 4.0, true);
        A.setValue(1, 1, 9.0, true);
        A.setValue(0, 1, 1.0, true);

        double x[2] = {2.0, 3.0};
        double y[2] = {0.0};

        A.gaxpy(x, y);

        double expected0 = 4.0 * 2.0 + 1.0 * 3.0;  // 8 + 3 = 11
        double expected1 = 1.0 * 2.0 + 9.0 * 3.0;  // 2 + 27 = 29

        REQUIRE(y[0] == Approx(expected0));
        REQUIRE(y[1] == Approx(expected1));
    }
}

TEST_CASE("CG Solver - Conjugate Gradient", "[solver]")
{
    SECTION("3x3 system solves correctly")
    {
        unsigned long n = 3;
        SymSparseProb prob(n);
        prob.setMaxiter(10000);
        prob.setToler(1e-15);

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
        prob.B[0] = 8.0;
        prob.B[1] = -2.0;
        prob.B[2] = 4.0;

        prob.CGSolve();

        // Verify solution: Ax = b
        double residual = 0.0;
        double b_norm = 0.0;
        for (unsigned long i = 0; i < n; i++) {
            b_norm += prob.B[i] * prob.B[i];
        }
        b_norm = std::sqrt(b_norm);

        // Compute Ax - b
        double Ax[3] = {0.0};
        prob.A.gaxpy(prob.X, Ax);
        for (unsigned long i = 0; i < n; i++) {
            double ri = Ax[i] - prob.B[i];
            residual += ri * ri;
        }
        residual = std::sqrt(residual);

        double relative_error = residual / b_norm;

        INFO("Solution: X = [" << prob.X[0] << ", " << prob.X[1] << ", " << prob.X[2] << "]");
        INFO("Relative residual: " << relative_error);

        REQUIRE(relative_error < 1e-10);
    }

    SECTION("JCG solver converges")
    {
        unsigned long n = 3;
        SymSparseProb prob(n);
        prob.setMaxiter(10000);
        prob.setToler(1e-15);

        prob.setMatrixValue(0, 0, 1.0, true);
        prob.setMatrixValue(1, 1, 2.0, true);
        prob.setMatrixValue(2, 2, 1.0, true);
        prob.setMatrixValue(0, 1, 5.0, true);
        prob.setMatrixValue(0, 2, -3.0, true);
        prob.setMatrixValue(1, 2, 2.0, true);

        prob.B[0] = 8.0;
        prob.B[1] = -2.0;
        prob.B[2] = 4.0;

        prob.JCGSolve();

        double residual = 0.0;
        double Ax[3] = {0.0};
        prob.A.gaxpy(prob.X, Ax);
        for (unsigned long i = 0; i < n; i++) {
            double ri = Ax[i] - prob.B[i];
            residual += ri * ri;
        }
        residual = std::sqrt(residual);

        REQUIRE(residual < 1e-8);

        REQUIRE(prob.X[0] == Approx(-2.0/3.0));
        REQUIRE(prob.X[1] == Approx(4.0/3.0));
        REQUIRE(prob.X[2] == Approx(-2.0/3.0));
    }
}

TEST_CASE("Matrix - Replace vs Add", "[matrix]")
{
    SECTION("replace=true overwrites value")
    {
        Sparse A(2);
        A.setValue(0, 0, 5.0, true);  // replace
        A.setValue(0, 0, 3.0, true);  // replace again

        REQUIRE(A.diag[0].val == 3.0);
    }

    SECTION("replace=false adds to value")
    {
        Sparse A(2);
        A.setValue(0, 0, 5.0, true);
        A.setValue(0, 0, 3.0, false);  // add

        REQUIRE(A.diag[0].val == 8.0);
    }
}
