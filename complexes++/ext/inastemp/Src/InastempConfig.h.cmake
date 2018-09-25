///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef INASTEMPCONFIG_H
#define INASTEMPCONFIG_H

// Define all macros (ADD-NEW-HERE)
#cmakedefine INASTEMP_USE_SCALAR

#cmakedefine INASTEMP_USE_SSE3

#cmakedefine INASTEMP_USE_SSSE3

#cmakedefine INASTEMP_USE_SSE41

#cmakedefine INASTEMP_USE_SSE42

#cmakedefine INASTEMP_USE_AVX

#cmakedefine INASTEMP_USE_AVX2

#cmakedefine INASTEMP_USE_AVX512COMMON

#cmakedefine INASTEMP_USE_AVX512KNL

#cmakedefine INASTEMP_USE_AVX512SKL

#cmakedefine INASTEMP_USE_ALTIVEC
#cmakedefine INASTEMP_USE_XL

// Inform about best one
#define INASTEMP_@INASTEMP_BESTTYPE@_IS_BEST
#define INASTEMP_BEST_TYPE @INASTEMP_BESTTYPE@

#ifndef INASTEMP_NO_BEST_INCLUDE
#include "@INASTEMP_BESTTYPE@/InaVec@INASTEMP_BESTTYPE@Float.hpp"
#include "@INASTEMP_BESTTYPE@/InaVec@INASTEMP_BESTTYPE@Double.hpp"
#else
template <class RealType>
class InaVec@INASTEMP_BESTTYPE@;
#endif

template <class RealType>
using InaVecBestType = InaVec@INASTEMP_BESTTYPE@<RealType>;

using InaVecBestTypeFloat = InaVec@INASTEMP_BESTTYPE@<float>;
using InaVecBestTypeDouble = InaVec@INASTEMP_BESTTYPE@<double>;


#endif
