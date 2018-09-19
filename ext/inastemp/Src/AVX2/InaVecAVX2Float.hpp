///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef INAVECAVX2FLOAT_HPP
#define INAVECAVX2FLOAT_HPP

#include "InastempConfig.h"
#include "AVX/InaVecAVXFloat.hpp"

#ifndef INASTEMP_USE_AVX2
#error InaVecAVX2<float> is included but AVX2 is not enable in the configuration
#endif

#include <tmmintrin.h>
#include <emmintrin.h>

template < class RealType >
class InaVecAVX2;

// AVX2 _mm256_abs_epi32 is not useful here
template <>
class alignas(32) InaVecAVX2< float > : public InaVecAVX< float > {
    using Parent = InaVecAVX< float >;

public:
    using InaVecAVX< float >::InaVecAVX;

    inline InaVecAVX2() {
    }

    inline InaVecAVX2(const InaVecAVX< float >& other)
    : Parent(other) {
    }

    inline static const char* GetName() {
        return "InaVecAVX2<float>";
    }

    inline static InaIfElse< InaVecAVX2< float > >::ThenClass If(const typename Parent::MaskType& inTest) {
        return InaIfElse< InaVecAVX2< float > >::IfClass().If(inTest);
    }
};


#endif
