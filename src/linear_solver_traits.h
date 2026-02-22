#ifndef LINEAR_SOLVER_TRAITS_H_INCLUDED
#define LINEAR_SOLVER_TRAITS_H_INCLUDED

#include <complex>
#include <cmath>

using namespace std;

struct RealTraits {
    using Scalar = double;
    
    static inline Scalar conjugate(const Scalar& x) { 
        return x; 
    }
    
    static inline Scalar zero() { 
        return 0.0; 
    }
    
    static inline double norm(const Scalar& x) { 
        return x * x; 
    }
    
    static inline double abs(const Scalar& x) { 
        return std::fabs(x); 
    }
    
    static constexpr bool is_complex = false;
    static constexpr bool supports_symmetric = true;
    static constexpr bool supports_hermitian = false;
    static constexpr bool enforces_symmetric_storage = true;
};

struct ComplexSymmetricTraits {
    using Scalar = std::complex<double>;
    
    static inline Scalar conjugate(const Scalar& x) { 
        return x;
    }
    
    static inline Scalar zero() { 
        return Scalar(0.0, 0.0); 
    }
    
    static inline double norm(const Scalar& x) { 
        return std::norm(x); 
    }
    
    static inline double abs(const Scalar& x) { 
        return std::abs(x); 
    }
    
    static constexpr bool is_complex = true;
    static constexpr bool supports_symmetric = true;
    static constexpr bool supports_hermitian = false;
    static constexpr bool enforces_symmetric_storage = true;
};

struct ComplexHermitianTraits {
    using Scalar = std::complex<double>;
    
    static inline Scalar conjugate(const Scalar& x) { 
        return std::conj(x);
    }
    
    static inline Scalar zero() { 
        return Scalar(0.0, 0.0); 
    }
    
    static inline double norm(const Scalar& x) { 
        return std::norm(x); 
    }
    
    static inline double abs(const Scalar& x) { 
        return std::abs(x); 
    }
    
    static constexpr bool is_complex = true;
    static constexpr bool supports_symmetric = false;
    static constexpr bool supports_hermitian = true;
    static constexpr bool enforces_symmetric_storage = false;
};

#endif
