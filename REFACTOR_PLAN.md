# Conjugate Gradients Solver Refactor Plan

## Overview

Refactor the codebase to use a **traits-based template architecture** that:
1. Eliminates ~70% code duplication between real and complex implementations
2. Provides compile-time validation of solver/matrix compatibility
3. Enables future extensions (Hermitian matrices, BiCGSTAB, GMRES)

---

## Current State

### Files (8 total, ~1200 lines)

```
src/
├── sparse_matrix.h           # Real sparse matrix (double)
├── sparse_matrix.cpp
├── sparse_matrix_complex.h   # Complex sparse matrix (complex<double>)
├── sparse_matrix_complex.cpp
├── linear_solver_real.h      # Real solvers (CG, JCG, SSOR-CG)
├── linear_solver_real.cpp
├── linear_solver_complex.h   # Complex solvers (COCG, COCG-SSOR)
└── linear_solver_complex.cpp
```

### Problems

| Issue | Description |
|-------|-------------|
| **Duplication** | ~95% identical code in sparse_matrix (real vs complex) |
| **Duplication** | ~80% identical code in linear_solver (CG vs COCG) |
| **No type safety** | Nothing prevents calling `solveWithCG()` on complex matrix |
| **Hard to extend** | Adding Hermitian support requires more duplicated files |
| **Memory bugs** | Manual `new`/`delete` scattered (some bugs found: `delete` vs `delete[]`) |

---

## Proposed Architecture

### 1. Traits Pattern

Traits encode scalar type properties at compile-time:

```cpp
// Real symmetric (SPD matrices)
struct RealTraits {
    using Scalar = double;
    static inline double conjugate(const double& x) { return x; }
    static inline double dot(const double& x, const double& y) { return x * y; }
    static inline double norm2(const double& x) { return x * x; }
    static constexpr bool is_complex = false;
    static constexpr bool supports_symmetric = true;
    static constexpr bool supports_hermitian = false;
};

// Complex symmetric (eddy current problems: A = A^T)
struct ComplexSymmetricTraits {
    using Scalar = std::complex<double>;
    static inline std::complex<double> conjugate(const std::complex<double>& x) { 
        return x; // NO conjugation - bilinear form x^T * y
    }
    static inline std::complex<double> dot(const std::complex<double>& x, 
                                           const std::complex<double>& y) { 
        return x * y; 
    }
    static inline double norm2(const std::complex<double>& x) { 
        return std::norm(x); 
    }
    static constexpr bool is_complex = true;
    static constexpr bool supports_symmetric = true;
    static constexpr bool supports_hermitian = false;
};

// Complex Hermitian (quantum mechanics: A = A^H) - FUTURE
struct ComplexHermitianTraits {
    using Scalar = std::complex<double>;
    static inline std::complex<double> conjugate(const std::complex<double>& x) { 
        return std::conj(x); // Conjugate - sesquilinear form x^H * y
    }
    static inline std::complex<double> dot(const std::complex<double>& x, 
                                           const std::complex<double>& y) { 
        return std::conj(x) * y; 
    }
    static inline double norm2(const std::complex<double>& x) { 
        return std::norm(x); 
    }
    static constexpr bool is_complex = true;
    static constexpr bool supports_symmetric = false;
    static constexpr bool supports_hermitian = true;
};
```

### 2. Solver/Matrix Compatibility

| Traits | Matrix Type | CG | COCG | BiCGSTAB | SSOR |
|--------|-------------|----|------|----------|------|
| `RealTraits` | Real SPD | ✓ | ✗ | ✓ | ✓ |
| `ComplexSymmetricTraits` | Complex Symmetric | ✗ | ✓ | ✓ | ✓ |
| `ComplexHermitianTraits` | Complex Hermitian | ✓ | ✗ | ✓ | ✓ |
| `RealGeneralTraits` (future) | General | ✗ | ✗ | ✓ | ✗ |

### 3. Templated Classes

```cpp
// Single sparse matrix for all scalar types
template<typename Traits>
class SparseMatrix {
    using Scalar = typename Traits::Scalar;
    
    struct Entry { 
        Scalar val; 
        unsigned long col; 
        Entry *next; 
    };
    
    struct DiagEntry { 
        Scalar val; 
        Entry *first; 
    };
    
    DiagEntry *diag;
    unsigned long N;
    
    void gaxpy(Scalar *x, Scalar *y);  // Works for all types
    void setValue(...);
};

// Single solver for all scalar types
template<typename Traits>
class LinearSolver {
    using Scalar = typename Traits::Scalar;
    
    SparseMatrix<Traits> matrix;
    Scalar *solution;
    Scalar *right_hand_side;
    
    // Compile-time validation via static_assert
    void solveWithCG() {
        static_assert(Traits::supports_symmetric || Traits::supports_hermitian,
                      "CG requires symmetric or Hermitian positive definite matrix");
        // ... implementation
    }
    
    void solveWithCOCG() {
        static_assert(Traits::is_complex && Traits::supports_symmetric,
                      "COCG requires complex symmetric matrix");
        // ... implementation
    }
    
    void solveWithBiCGSTAB() {
        // Available for all types (future implementation)
    }
    
    void solveWithSSOR(double omega);
    void solveWithJCG();  // Real only (static_assert)
};
```

### 4. Type Aliases (Backward Compatibility)

```cpp
// Keep existing names for backward compatibility
using Sparse = SparseMatrix<RealTraits>;
using ComplexSparseMatrix = SparseMatrix<ComplexSymmetricTraits>;

using RealLinearSolver = LinearSolver<RealTraits>;
using ComplexLinearSolver = LinearSolver<ComplexSymmetricTraits>;

// New names for clarity (optional)
using RealSymmetricSolver = RealLinearSolver;
using ComplexSymmetricSolver = ComplexLinearSolver;
using ComplexHermitianSolver = LinearSolver<ComplexHermitianTraits>;
```

### 5. Usage Examples

```cpp
// Real SPD problem (existing code works unchanged)
RealLinearSolver solver(100);
solver.setMatrixValue(0, 0, 1.0, true);
solver.solveWithCG();      // ✓ compiles
solver.solveWithCOCG();    // ✗ compile error: "COCG requires complex symmetric"

// Complex symmetric (eddy current - existing code works unchanged)
ComplexLinearSolver solver(100);
solver.setMatrixValue(0, 0, {1.0, 0.0}, true);
solver.solveWithCOCG();    // ✓ compiles
solver.solveWithCG();      // ✗ compile error: "CG requires symmetric or Hermitian"

// Complex Hermitian (NEW - future use case)
ComplexHermitianSolver solver(100);
solver.solveWithCG();      // ✓ compiles (CG works for Hermitian PD)
solver.solveWithBiCGSTAB(); // ✓ also works

// General non-symmetric (FUTURE)
using RealGeneralSolver = LinearSolver<RealGeneralTraits>;
RealGeneralSolver solver(100);
solver.solveWithBiCGSTAB(); // ✓ works for any matrix
solver.solveWithCG();       // ✗ compile error
```

---

## Implementation Phases

### Phase 1: Create Traits Header

**File:** `src/linear_solver_traits.h`

- `RealTraits` struct
- `ComplexSymmetricTraits` struct
- `ComplexHermitianTraits` struct (forward declaration, not fully used yet)
- `SolverType` enum for runtime debugging

### Phase 2: Create Templated Sparse Matrix

**File:** `src/sparse_matrix.h`

- `template<typename Traits> class SparseMatrix`
- Internal `Entry` and `DiagEntry` structs (templated)
- `setValue(...)` - replaces value (default behavior)
- `addToValue(...)` - adds to existing value (clearer than `bool replace`)
- `gaxpy()`, `printOut()` methods
- Type aliases at bottom for backward compatibility

### Phase 2b: API Change - Remove `bool replace` [DONE]

**Completed:** Split into two clear methods:

**Before:**
```cpp
void setValue(unsigned long i, unsigned long j, double valu, bool replace);
```

**After:**
```cpp
void setValue(unsigned long i, unsigned long j, Scalar valu);    // Replace
void addToValue(unsigned long i, unsigned long j, Scalar valu);  // Add
```

**Migration:**
- Solver code: `setMatrixValue(..., replace=true)` → `setMatrixValue()` (unchanged)
- Test code: `setValue(0, 0, 3.0, false)` → `addToValue(0, 0, 3.0)`

**Files modified:**
- `src/sparse_matrix.h` - Added `addToValue()`, removed `bool replace`
- `src/sparse_matrix.cpp` - Split implementation
- `src/sparse_matrix_complex.h` - Added `addToValue()`, removed `bool replace`
- `src/sparse_matrix_complex.cpp` - Split implementation
- `src/linear_solver_real.h` - Added `addToMatrixValue()`
- `src/linear_solver_real.cpp` - Added `addToMatrixValue()` impl
- `src/linear_solver_complex.h` - Added `addToMatrixValue()`
- `src/linear_solver_complex.cpp` - Added `addToMatrixValue()` impl
- `tests/test_solvers.cpp` - Updated to use new API

### Phase 3: Create Templated Linear Solver

**File:** `src/linear_solver.h`

- `template<typename Traits> class LinearSolver`
- `solveWithCG()` - with `static_assert` for Real/Hermitian
- `solveWithCOCG()` - with `static_assert` for ComplexSymmetric
- `solveWithSSOR()` - preconditioned CG/COCG
- `solveWithJCG()` - Jacobi preconditioned (real only)
- `solveWithBiCGSTAB()` - stub for future
- `getIterations()`, `setTolerance()`, etc.

### Phase 4: Update Build System

**File:** `src/CMakeLists.txt`

```cmake
# Before
set(SOURCES
    sparse_matrix.cpp
    sparse_matrix_complex.cpp
    linear_solver_real.cpp
    linear_solver_complex.cpp
)

# After
set(SOURCES
    # No .cpp files - templates are header-only
)
```

### Phase 5: Delete Old Files

```bash
# Delete old implementations
rm src/sparse_matrix.h
rm src/sparse_matrix.cpp
rm src/sparse_matrix_complex.h
rm src/sparse_matrix_complex.cpp
rm src/linear_solver_real.h
rm src/linear_solver_real.cpp
rm src/linear_solver_complex.h
rm src/linear_solver_complex.cpp
```

### Phase 6: Update Tests

**File:** `tests/test_solvers.cpp`

```cpp
// Change includes
#include "linear_solver_real.h"      → #include "linear_solver.h"
#include "linear_solver_complex.h"   → #include "linear_solver.h"

// No other changes needed - type aliases preserve API
```

---

## File Structure After Refactor

```
src/
├── CMakeLists.txt              # Modified
├── linear_solver_traits.h      # NEW - traits definitions
├── sparse_matrix.h             # NEW - templated matrix
├── linear_solver.h             # NEW - templated solver
└── (old files deleted)

tests/
└── test_solvers.cpp            # Modified includes only
```

**Before:** 8 files, ~1200 lines  
**After:** 3 files, ~500 lines  
**Reduction:** ~58% code reduction

---

## Benefits

| Benefit | Description |
|---------|-------------|
| **Code reduction** | 8 files → 3 files, ~1200 lines → ~500 lines |
| **Type safety** | Wrong solver for matrix type = compile error |
| **Single source of truth** | Fix bug once, applies to all scalar types |
| **Extensibility** | Add `FloatTraits`, `QuadTraits`, `GeneralTraits` easily |
| **Performance** | Zero runtime overhead (all compile-time) |
| **Maintainability** | No more sync issues between real/complex versions |
| **Future-proof** | Hermitian, BiCGSTAB, GMRES add once, work for all |

---

## Migration Risks & Mitigation

| Risk | Mitigation |
|------|------------|
| Template compile errors | Keep type aliases for backward compatibility |
| Tests fail | Run existing 5 tests, verify all pass |
| IntelliSense breaks | `compile_commands.json` already configured |
| Performance regression | Templates optimize same as concrete types |

---

## Future Extensions

### Hermitian Support (Phase 7)
```cpp
ComplexHermitianSolver solver(n);
solver.solveWithCG();  // Uses conjugate dot product
```

### BiCGSTAB (Phase 8)
```cpp
template<typename Traits>
void LinearSolver<Traits>::solveWithBiCGSTAB() {
    // Single implementation works for ALL matrix types
    // No symmetry required
}
```

### GMRES (Phase 9)
```cpp
template<typename Traits>
void LinearSolver<Traits>::solveWithGMRES(unsigned long restart) {
    // Krylov subspace method for general matrices
}
```

### Block Solvers (Phase 10)
```cpp
template<typename Traits, unsigned long BlockSize>
class BlockLinearSolver {
    // Solve multiple RHS simultaneously
};
```

---

## References

- [Traits Pattern](https://en.wikibooks.org/wiki/More_C%2B%2B_Idioms/Traits)
- [SFINAE and static_assert](https://en.cppreference.com/w/cpp/language/static_assert)
- [COCG Algorithm](https://epubs.siam.org/doi/10.1137/0914053)
- [BiCGSTAB Paper](https://epubs.siam.org/doi/10.1137/0914041)
