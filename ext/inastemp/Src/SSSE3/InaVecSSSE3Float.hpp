///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef INAVECSSSE3FLOAT_HPP
#define INAVECSSSE3FLOAT_HPP

#include "InastempConfig.h"
#include "SSE3/InaVecSSE3Float.hpp"

#ifndef INASTEMP_USE_SSSE3
#error InaVecSSSE3<float> is included but SSSE3 is not enable in the configuration
#endif

#include <tmmintrin.h>
#include <emmintrin.h>

// Forward declarations
template < class RealType >
class InaVecSSSE3;

// SSSE3 add _mm_sign_epi32/_mm_abs_epi32 but not useful here
template <>
class alignas(16) InaVecSSSE3< float > : public InaVecSSE3< float > {
    using Parent = InaVecSSE3< float >;

public:
    using InaVecSSE3< float >::InaVecSSE3;

    inline InaVecSSSE3() {
    }

    inline InaVecSSSE3(const InaVecSSE3< float >& other)
    : Parent(other) {
    }

    inline static const char* GetName() {
        return "InaVecSSSE3<float>";
    }

    inline static InaIfElse< InaVecSSSE3< float > >::ThenClass If(const typename Parent::MaskType& inTest) {
        return InaIfElse< InaVecSSSE3< float > >::IfClass().If(inTest);
    }
};

#endif
