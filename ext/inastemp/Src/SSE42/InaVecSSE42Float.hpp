///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef INAVECSSE42FLOAT_HPP
#define INAVECSSE42FLOAT_HPP

#include "InastempConfig.h"
#include "SSE41/InaVecSSE41Float.hpp"

#ifndef INASTEMP_USE_SSE42
#error InaVecSSE42<float> is included but SSE42 is not enable in the configuration
#endif

#include <tmmintrin.h>
#include <emmintrin.h>

template < class RealType >
class InaVecSSE42;

template <>
class alignas(16) InaVecSSE42< float > : public InaVecSSE41< float > {
    using Parent = InaVecSSE41< float >;

public:
    using InaVecSSE41< float >::InaVecSSE41;

    inline InaVecSSE42() {
    }

    inline InaVecSSE42(const InaVecSSE41< float >& other)
    : Parent(other) {
    }

    inline static const char* GetName() {
        return "InaVecSSE42<float>";
    }

    inline static InaIfElse< InaVecSSE42< float > >::ThenClass If(const typename Parent::MaskType& inTest) {
        return InaIfElse< InaVecSSE42< float > >::IfClass().If(inTest);
    }
};


#endif
