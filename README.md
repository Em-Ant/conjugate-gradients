# Conjugate Gradients

### Simple Implementation of the Iterative Solver

See [this paper](https://www.cs.cmu.edu/~quake-papers/painless-conjugate-gradient.pdf) for a great explanation of the subject.

It's suitable for Symmetric Linear (Sparse) Problems.

## Architecture

- **Traits-based template design** eliminates 70% code duplication
- **Compile-time validation** prevents wrong solver/matrix combinations
- **Backward compatible** with existing API
- **Extensible** for Hermitian, BiCGSTAB, GMRES

## Features

- **Real solvers**: CG, JCG, SSOR-CG
- **Complex solvers**: COCG, COCG-SSOR (eddy current problems)
- **Preconditioners**: Jacobi, SSOR
- **Matrix storage**: Symmetric List Of Lists (i ≤ j enforced)
- **Type safety**: Wrong solver = compile error
- **Performance**: Zero runtime overhead (all compile-time)

## Usage

```cpp
// Real SPD problem
RealLinearSolver solver(100);
solver.setMatrixValue(0, 0, 1.0);
solver.solveWithCG();      // ✓ compiles
solver.solveWithCOCG();    // ✗ compile error

// Complex symmetric (eddy current)
ComplexLinearSolver solver(100);
solver.setMatrixValue(0, 0, {1.0, 0.0});
solver.solveWithCOCG();    // ✓ compiles
solver.solveWithCG();      // ✗ compile error
```

## Build

```bash
mkdir build && cd build
cmake ..
make
```

## Test

```bash
./tests/test_solvers
```

## Mathematical Background

The system matrix is modeled as a **List Of Lists**.

It features two simple preconditioners:

- [Jacobi](https://en.wikipedia.org/wiki/Preconditioner#Jacobi_.28or_diagonal.29_preconditioner)
- [SSOR](https://en.wikipedia.org/wiki/Successive_over-relaxation)

I used it as a solver for simple 2D [FEM](https://en.wikipedia.org/wiki/Finite_element_method) [Poisson Problems](https://en.wikipedia.org/wiki/Poisson's_equation).
