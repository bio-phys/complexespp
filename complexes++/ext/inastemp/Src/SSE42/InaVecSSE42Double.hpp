///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef INAVECSSE42DOUBLE_HPP
#define INAVECSSE42DOUBLE_HPP

#include "InastempConfig.h"
#include "SSE41/InaVecSSE41Double.hpp"

#ifndef INASTEMP_USE_SSE42
#error InaVecSSE42<double> is included but SSE42 is not enable in the configuration
#endif

#include <tmmintrin.h>
#include <emmintrin.h>

template < class RealType >
class InaVecSSE42;

template <>
class alignas(16) InaVecSSE42< double > : public InaVecSSE41< double > {
    using Parent = InaVecSSE41< double >;

public:
    using InaVecSSE41< double >::InaVecSSE41;

    inline InaVecSSE42() {
    }

    inline InaVecSSE42(const InaVecSSE41< double >& other)
    : Parent(other) {
    }

    inline static const char* GetName() {
        return "InaVecSSE42<double>";
    }

    inline static InaIfElse< InaVecSSE42< double > >::ThenClass If(const typename Parent::MaskType& inTest) {
        return InaIfElse< InaVecSSE42< double > >::IfClass().If(inTest);
    }
};

#endif
