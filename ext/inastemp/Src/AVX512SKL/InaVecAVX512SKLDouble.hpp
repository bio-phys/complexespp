///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef INAVECAVX512SKLDOUBLE_HPP
#define INAVECAVX512SKLDOUBLE_HPP

#include "InastempConfig.h"
#include "Common/InaIfElse.hpp"
#include "Common/InaUtils.hpp"
#include "AVX512COMMON/InaVecAVX512COMMONDouble.hpp"

#ifndef INASTEMP_USE_AVX512SKL
#error InaVecAVX512SKL<double> is included but AVX512SKL is not enable in the configuration
#endif

#include "Common/InaFastExp.hpp"

#include <immintrin.h>

#include <cmath>

// Forward declarations
template < class RealType >
using InaVecMaskAVX512SKL = InaVecMaskAVX512COMMON< RealType >;

template < class RealType >
using InaVecAVX512SKL = InaVecAVX512COMMON< RealType >;


#endif
