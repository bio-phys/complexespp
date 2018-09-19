///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef INAVECSSSE3DOUBLE_HPP
#define INAVECSSSE3DOUBLE_HPP

#include "SSE3/InaVecSSE3Double.hpp"

#ifndef INASTEMP_USE_SSSE3
#error InaVecSSSE3<double> is included but SSSE3 is not enable in the configuration
#endif

#include <tmmintrin.h>
#include <emmintrin.h>

// Forward declarations
template < class RealType >
class InaVecSSSE3;

// No optimization to make for SSSE in double
template <>
class alignas(16) InaVecSSSE3< double > : public InaVecSSE3< double > {
    using Parent = InaVecSSE3< double >;

public:
    using InaVecSSE3< double >::InaVecSSE3;

    inline InaVecSSSE3() {
    }

    inline InaVecSSSE3(const InaVecSSE3< double >& other)
    : Parent(other) {
    }

    inline static const char* GetName() {
        return "InaVecSSSE3<double>";
    }

    inline static InaIfElse< InaVecSSSE3< double > >::ThenClass If(const typename Parent::MaskType& inTest) {
        return InaIfElse< InaVecSSSE3< double > >::IfClass().If(inTest);
    }
};

#endif
