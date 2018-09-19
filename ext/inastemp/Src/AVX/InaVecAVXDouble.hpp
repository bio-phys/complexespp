///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef INAVECAVXDOUBLE_HPP
#define INAVECAVXDOUBLE_HPP

#include "InastempConfig.h"
#include "Common/InaIfElse.hpp"
#include "Common/InaUtils.hpp"

#ifndef INASTEMP_USE_AVX
#error InaVecAVX<double> is included but AVX is not enable in the configuration
#endif

#include "Common/InaFastExp.hpp"

#include <immintrin.h>
#include <cmath>
#include <initializer_list>

// Forward declarations
template < class RealType >
class InaVecMaskAVX;

template < class RealType >
class InaVecAVX;

// Mask type
template <>
class alignas(32) InaVecMaskAVX< double > {
    __m256i mask;

public:
    // Classic constructors
    inline InaVecMaskAVX() {
    }

    inline InaVecMaskAVX(const InaVecMaskAVX&) = default;
    inline InaVecMaskAVX& operator=(const InaVecMaskAVX&) = default;

    // Native data type compatibility
    inline /*not explicit*/ InaVecMaskAVX(const __m256i inMask)
    : mask(inMask) {
    }

    inline InaVecMaskAVX& operator=(const __m256i inMask) {
        mask = inMask;
        return (*this);
    }

    inline explicit operator __m256i() const {
        return mask;
    }

    inline __m256i getMask() const {
        return mask;
    }

    // Bool data type compatibility
    inline explicit InaVecMaskAVX(const bool inBool) {
        mask = (inBool ? _mm256_set1_epi64x(static_cast< long long >(0xFFFFFFFFFFFFFFFFL)) : _mm256_setzero_si256());
    }

    inline InaVecMaskAVX& operator=(const bool inBool) {
        mask = (inBool ? _mm256_set1_epi64x(static_cast< long long >(0xFFFFFFFFFFFFFFFFL)) : _mm256_setzero_si256());
        return (*this);
    }

    // Binary methods
    inline InaVecMaskAVX Not() const {
        return NotAnd(mask, _mm256_set1_epi64x(static_cast< long long >(0xFFFFFFFFFFFFFFFFL)));
    }

    inline bool isAllTrue() const {
        // true if all FF => !FF => 0 & FF => 0
        return _mm256_testc_si256(mask, _mm256_set1_epi64x(static_cast< long long >(0xFFFFFFFFFFFFFFFFL)));
    }

    inline bool isAllFalse() const {
        // true if all zero
        return _mm256_testz_si256(mask, mask);
    }

    // Double args methods
    inline static InaVecMaskAVX And(const InaVecMaskAVX& inMask1, const InaVecMaskAVX& inMask2) {
        // AVX2 return InaVecMaskAVX(_mm256_and_si256(inMask1.mask, inMask2.mask));
        return InaVecMaskAVX(_mm256_castpd_si256(_mm256_and_pd(_mm256_castsi256_pd(inMask1.mask),
                                                               _mm256_castsi256_pd(inMask2.mask))));
    }

    inline static InaVecMaskAVX NotAnd(const InaVecMaskAVX& inMask1, const InaVecMaskAVX& inMask2) {
        // AVX2 return InaVecMaskAVX(_mm256_andnot_si256(inMask1.mask, inMask2.mask));
        return InaVecMaskAVX(_mm256_castpd_si256(_mm256_andnot_pd(_mm256_castsi256_pd(inMask1.mask),
                                                                  _mm256_castsi256_pd(inMask2.mask))));
    }

    inline static InaVecMaskAVX Or(const InaVecMaskAVX& inMask1, const InaVecMaskAVX& inMask2) {
        // AVX2 return InaVecMaskAVX(_mm256_or_si256(inMask1.mask, inMask2.mask));
        return InaVecMaskAVX(_mm256_castpd_si256(_mm256_or_pd(_mm256_castsi256_pd(inMask1.mask),
                                                              _mm256_castsi256_pd(inMask2.mask))));
    }

    inline static InaVecMaskAVX Xor(const InaVecMaskAVX& inMask1, const InaVecMaskAVX& inMask2) {
        // AVX2 return InaVecMaskAVX(_mm256_xor_si256(inMask1.mask, inMask2.mask));
        return InaVecMaskAVX(_mm256_castpd_si256(_mm256_xor_pd(_mm256_castsi256_pd(inMask1.mask),
                                                               _mm256_castsi256_pd(inMask2.mask))));
    }

    inline static bool IsEqual(const InaVecMaskAVX& inMask1, const InaVecMaskAVX& inMask2) {
        return _mm256_testz_si256(_mm256_castpd_si256(_mm256_xor_pd(_mm256_castsi256_pd(inMask1.mask),
                                                                    _mm256_castsi256_pd(inMask2.mask))),
                                  _mm256_set1_epi64x(static_cast< long long >(0xFFFFFFFFFFFFFFFFL))); // return CF
    }

    inline static bool IsNotEqual(const InaVecMaskAVX& inMask1, const InaVecMaskAVX& inMask2) {
        return !_mm256_testz_si256(_mm256_castpd_si256(_mm256_xor_pd(_mm256_castsi256_pd(inMask1.mask),
                                                                     _mm256_castsi256_pd(inMask2.mask))),
                                   _mm256_set1_epi64x(static_cast< long long >(0xFFFFFFFFFFFFFFFFL))); // return CF
    }
};

// Mask must have operators
inline InaVecMaskAVX< double > operator&(const InaVecMaskAVX< double >& inMask1, const InaVecMaskAVX< double >& inMask2) {
    return InaVecMaskAVX< double >::And(inMask1, inMask2);
}

inline InaVecMaskAVX< double > operator|(const InaVecMaskAVX< double >& inMask1, const InaVecMaskAVX< double >& inMask2) {
    return InaVecMaskAVX< double >::Or(inMask1, inMask2);
}

inline InaVecMaskAVX< double > operator^(const InaVecMaskAVX< double >& inMask1, const InaVecMaskAVX< double >& inMask2) {
    return InaVecMaskAVX< double >::Xor(inMask1, inMask2);
}

inline bool operator==(const InaVecMaskAVX< double >& inMask1, const InaVecMaskAVX< double >& inMask2) {
    return InaVecMaskAVX< double >::IsEqual(inMask1, inMask2);
}

inline bool operator!=(const InaVecMaskAVX< double >& inMask1, const InaVecMaskAVX< double >& inMask2) {
    return InaVecMaskAVX< double >::IsNotEqual(inMask1, inMask2);
}

// Vec type
template <>
class alignas(32) InaVecAVX< double > {
protected:
    __m256d vec;

public:
    using VecRawType            = __m256d;
    using MaskType              = InaVecMaskAVX< double >;
    using RealType              = double;
    static const int VecLength  = 4;
    static const int Alignement = 32;

    inline InaVecAVX() {
    }
    inline InaVecAVX(const InaVecAVX&) = default;
    inline InaVecAVX& operator=(const InaVecAVX&) = default;

    // Constructor from raw type
    inline /*not explicit*/ InaVecAVX(const __m256d inVec)
    : vec(inVec) {
    }

    inline InaVecAVX& operator=(const __m256d inVec) {
        vec = inVec;
        return *this;
    }

    inline void setFromRawType(const __m256d inVec) {
        vec = inVec;
    }

    inline explicit operator __m256d() const {
        return vec;
    }

    inline __m256d getVec() const {
        return vec;
    }

    // Constructor from scalar
    inline /*not explicit*/ InaVecAVX(const double val)
    : vec(_mm256_set1_pd(val)) {
    }

    inline InaVecAVX& operator=(const double val) {
        vec = _mm256_set1_pd(val);
        return *this;
    }

    inline void setFromScalar(const double val) {
        vec = _mm256_set1_pd(val);
    }

    // Constructor from vec
    inline InaVecAVX(const std::initializer_list< double > lst)
    : InaVecAVX(lst.begin()) {
    }

    inline explicit InaVecAVX(const double ptr[])
    : vec(_mm256_loadu_pd(ptr)) {
    }

    inline InaVecAVX& setFromArray(const double ptr[]) {
        vec = _mm256_loadu_pd(ptr);
        return *this;
    }

    inline InaVecAVX& setFromAlignedArray(const double ptr[]) {
        vec = _mm256_load_pd(ptr);
        return *this;
    }

    inline InaVecAVX& setFromIndirectArray(const double values[], const int inIndirection[]) {
        vec = _mm256_set_pd(
            values[inIndirection[3]],
            values[inIndirection[2]],
            values[inIndirection[1]],
            values[inIndirection[0]]);
        return *this;
    }

    inline InaVecAVX& setFromIndirect2DArray(const double inArray[], const int inIndirection1[],
                                             const int inLeadingDimension, const int inIndirection2[]) {
        vec = _mm256_set_pd(
            inArray[inIndirection1[3] * inLeadingDimension + inIndirection2[3]],
            inArray[inIndirection1[2] * inLeadingDimension + inIndirection2[2]],
            inArray[inIndirection1[1] * inLeadingDimension + inIndirection2[1]],
            inArray[inIndirection1[0] * inLeadingDimension + inIndirection2[0]]);
        return *this;
    }

    // Move back to array
    inline void storeInArray(double ptr[]) const {
        _mm256_storeu_pd(ptr, vec);
    }

    inline void storeInAlignedArray(double ptr[]) const {
        _mm256_store_pd(ptr, vec);
    }

    // Acce to individual values
    inline double at(const int index) const {
        alignas(Alignement) double allval[VecLength];
        _mm256_store_pd(allval, vec);
        return allval[index];
    }

    // Horizontal operation
    inline double horizontalSum() const {
        const __m128d valupper = _mm256_extractf128_pd(vec, 1);
        const __m128d rest     = _mm256_castpd256_pd128(vec);
        const __m128d valval   = _mm_add_pd(valupper, rest);
        const __m128d res      = _mm_add_pd(_mm_permute_pd(valval, 1), valval);
        return _mm_cvtsd_f64(res);
    }

    inline double horizontalMul() const {
        const __m128d valupper = _mm256_extractf128_pd(vec, 1);
        const __m128d rest     = _mm256_castpd256_pd128(vec);
        const __m128d valval   = _mm_mul_pd(valupper, rest);
        const __m128d res      = _mm_mul_pd(_mm_permute_pd(valval, 1), valval);
        return _mm_cvtsd_f64(res);
    }

    inline InaVecAVX sqrt() const {
        return _mm256_sqrt_pd(vec);
    }

    inline InaVecAVX exp() const {
#ifdef __INTEL_COMPILER
        return _mm256_exp_pd(vec);
#else
        const __m256d COEFF_LOG2E = _mm256_set1_pd(double(InaFastExp::CoeffLog2E()));
        const __m256d COEFF_A     = _mm256_set1_pd(double(InaFastExp::CoeffA64()));
        const __m256d COEFF_B     = _mm256_set1_pd(double(InaFastExp::CoeffB64()));
        const __m256d COEFF_P5_X  = _mm256_set1_pd(double(InaFastExp::GetCoefficient9_8()));
        const __m256d COEFF_P5_Y  = _mm256_set1_pd(double(InaFastExp::GetCoefficient9_7()));
        const __m256d COEFF_P5_Z  = _mm256_set1_pd(double(InaFastExp::GetCoefficient9_6()));
        const __m256d COEFF_P5_A  = _mm256_set1_pd(double(InaFastExp::GetCoefficient9_5()));
        const __m256d COEFF_P5_B  = _mm256_set1_pd(double(InaFastExp::GetCoefficient9_4()));
        const __m256d COEFF_P5_C  = _mm256_set1_pd(double(InaFastExp::GetCoefficient9_3()));
        const __m256d COEFF_P5_D  = _mm256_set1_pd(double(InaFastExp::GetCoefficient9_2()));
        const __m256d COEFF_P5_E  = _mm256_set1_pd(double(InaFastExp::GetCoefficient9_1()));
        const __m256d COEFF_P5_F  = _mm256_set1_pd(double(InaFastExp::GetCoefficient9_0()));

        __m256d x = _mm256_mul_pd(vec, COEFF_LOG2E);

        const __m256d fractional_part = _mm256_sub_pd(x, InaVecAVX(x).floor().vec);

        __m256d factor = _mm256_add_pd(_mm256_mul_pd(_mm256_add_pd(
                                                         _mm256_mul_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_add_pd(
                                                                                                       _mm256_mul_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_add_pd(
                                                                                                                                                     _mm256_mul_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_add_pd(_mm256_mul_pd(
                                                                                                                                                                                                                 COEFF_P5_X, fractional_part),
                                                                                                                                                                                                             COEFF_P5_Y),
                                                                                                                                                                                               fractional_part),
                                                                                                                                                                                 COEFF_P5_Z),
                                                                                                                                                                   fractional_part),
                                                                                                                                                     COEFF_P5_A),
                                                                                                                                                 fractional_part),
                                                                                                                                   COEFF_P5_B),
                                                                                                                     fractional_part),
                                                                                                       COEFF_P5_C),
                                                                                                   fractional_part),
                                                                                     COEFF_P5_D),
                                                                       fractional_part),
                                                         COEFF_P5_E),
                                                     fractional_part),
                                       COEFF_P5_F);

        x = _mm256_sub_pd(x, factor);

        x = _mm256_add_pd(_mm256_mul_pd(COEFF_A, x), COEFF_B);

        __m128d valupper = _mm256_extractf128_pd(x, 1);
        __m128d vallower = _mm256_castpd256_pd128(x);

        alignas(64) long long int allvalint[VecLength] = { _mm_cvtsd_si64(vallower),
                                                           _mm_cvtsd_si64(_mm_shuffle_pd(vallower, vallower, 1)),
                                                           _mm_cvtsd_si64(valupper),
                                                           _mm_cvtsd_si64(_mm_shuffle_pd(valupper, valupper, 1)) };

        return _mm256_castsi256_pd(_mm256_load_si256(reinterpret_cast< const __m256i* >(allvalint)));
#endif
    }

    inline InaVecAVX expLowAcc() const {
        const __m256d COEFF_LOG2E = _mm256_set1_pd(double(InaFastExp::CoeffLog2E()));
        const __m256d COEFF_A     = _mm256_set1_pd(double(InaFastExp::CoeffA64()));
        const __m256d COEFF_B     = _mm256_set1_pd(double(InaFastExp::CoeffB64()));
        const __m256d COEFF_P5_C  = _mm256_set1_pd(double(InaFastExp::GetCoefficient4_3()));
        const __m256d COEFF_P5_D  = _mm256_set1_pd(double(InaFastExp::GetCoefficient4_2()));
        const __m256d COEFF_P5_E  = _mm256_set1_pd(double(InaFastExp::GetCoefficient4_1()));
        const __m256d COEFF_P5_F  = _mm256_set1_pd(double(InaFastExp::GetCoefficient4_0()));

        __m256d x = _mm256_mul_pd(vec, COEFF_LOG2E);

        const __m256d fractional_part = _mm256_sub_pd(x, InaVecAVX(x).floor().vec);

        __m256d factor = _mm256_add_pd(_mm256_mul_pd(_mm256_add_pd(
                                                         _mm256_mul_pd(_mm256_add_pd(_mm256_mul_pd(
                                                                                         COEFF_P5_C, fractional_part),
                                                                                     COEFF_P5_D),
                                                                       fractional_part),
                                                         COEFF_P5_E),
                                                     fractional_part),
                                       COEFF_P5_F);

        x = _mm256_sub_pd(x, factor);

        x = _mm256_add_pd(_mm256_mul_pd(COEFF_A, x), COEFF_B);

        __m128d valupper = _mm256_extractf128_pd(x, 1);
        __m128d vallower = _mm256_castpd256_pd128(x);

        alignas(64) long long int allvalint[VecLength] = { _mm_cvtsd_si64(vallower),
                                                           _mm_cvtsd_si64(_mm_shuffle_pd(vallower, vallower, 1)),
                                                           _mm_cvtsd_si64(valupper),
                                                           _mm_cvtsd_si64(_mm_shuffle_pd(valupper, valupper, 1)) };

        return _mm256_castsi256_pd(_mm256_load_si256(reinterpret_cast< const __m256i* >(allvalint)));
    }

    inline InaVecAVX rsqrt() const {
        return _mm256_set1_pd(1) / _mm256_sqrt_pd(vec);
    }

    inline InaVecAVX abs() const {
        const __m256d minus0 = _mm256_castsi256_pd(_mm256_set1_epi64x(static_cast< long long >(0x8000000000000000L)));
        return _mm256_andnot_pd(minus0, vec);
    }

    inline InaVecAVX floor() const {
        return _mm256_floor_pd(vec);
    }

    inline InaVecAVX signOf() const {
        const __m256d minus0 = _mm256_castsi256_pd(_mm256_set1_epi64x(static_cast< long long >(0x8000000000000000L)));
        const __m256d signs  = _mm256_and_pd(vec, minus0);
        return _mm256_and_pd(_mm256_cmp_pd(_mm256_setzero_pd(), vec, _CMP_NEQ_OQ),
                             _mm256_or_pd(signs, _mm256_set1_pd(1)));
    }

    inline InaVecAVX isPositive() const {
        const __m256d greater = _mm256_cmp_pd(_mm256_setzero_pd(), vec, _CMP_LE_OQ);
        const __m256d ones    = _mm256_set1_pd(1);
        return _mm256_and_pd(greater, ones);
    }

    inline InaVecAVX isNegative() const {
        const __m256d less = _mm256_cmp_pd(_mm256_setzero_pd(), vec, _CMP_GE_OQ);
        const __m256d ones = _mm256_set1_pd(1);
        return _mm256_and_pd(less, ones);
    }

    inline InaVecAVX isPositiveStrict() const {
        const __m256d greater = _mm256_cmp_pd(_mm256_setzero_pd(), vec, _CMP_LT_OQ);
        const __m256d ones    = _mm256_set1_pd(1);
        return _mm256_and_pd(greater, ones);
    }

    inline InaVecAVX isNegativeStrict() const {
        const __m256d less = _mm256_cmp_pd(_mm256_setzero_pd(), vec, _CMP_GT_OQ);
        const __m256d ones = _mm256_set1_pd(1);
        return _mm256_and_pd(less, ones);
    }

    inline InaVecAVX isZero() const {
        const __m256d equalZero = _mm256_cmp_pd(_mm256_setzero_pd(), vec, _CMP_EQ_OQ);
        const __m256d ones      = _mm256_set1_pd(1);
        return _mm256_and_pd(equalZero, ones);
    }

    inline InaVecAVX isNotZero() const {
        const __m256d equalZero = _mm256_cmp_pd(_mm256_setzero_pd(), vec, _CMP_NEQ_OQ);
        const __m256d ones      = _mm256_set1_pd(1);
        return _mm256_and_pd(equalZero, ones);
    }

    inline InaVecMaskAVX< double > isPositiveMask() const {
        return _mm256_castpd_si256(_mm256_cmp_pd(_mm256_setzero_pd(), vec, _CMP_LE_OQ));
    }

    inline InaVecMaskAVX< double > isNegativeMask() const {
        return _mm256_castpd_si256(_mm256_cmp_pd(_mm256_setzero_pd(), vec, _CMP_GE_OQ));
    }

    inline InaVecMaskAVX< double > isPositiveStrictMask() const {
        return _mm256_castpd_si256(_mm256_cmp_pd(_mm256_setzero_pd(), vec, _CMP_LT_OQ));
    }

    inline InaVecMaskAVX< double > isNegativeStrictMask() const {
        return _mm256_castpd_si256(_mm256_cmp_pd(_mm256_setzero_pd(), vec, _CMP_GT_OQ));
    }

    inline InaVecMaskAVX< double > isZeroMask() const {
        return _mm256_castpd_si256(_mm256_cmp_pd(_mm256_setzero_pd(), vec, _CMP_EQ_OQ));
    }

    inline InaVecMaskAVX< double > isNotZeroMask() const {
        return _mm256_castpd_si256(_mm256_cmp_pd(_mm256_setzero_pd(), vec, _CMP_NEQ_OQ));
    }

    // Static basic methods
    inline static InaVecAVX GetZero() {
        return InaVecAVX(_mm256_setzero_pd());
    }

    inline static InaVecAVX GetOne() {
        return InaVecAVX(_mm256_set1_pd(1));
    }

    inline static InaVecAVX Min(const InaVecAVX& inVec1, const InaVecAVX& inVec2) {
        return _mm256_min_pd(inVec1.vec, inVec2.vec);
    }

    inline static InaVecAVX Max(const InaVecAVX& inVec1, const InaVecAVX& inVec2) {
        return _mm256_max_pd(inVec1.vec, inVec2.vec);
    }

    inline static InaVecAVX IsLowerOrEqual(const InaVecAVX& inVec1, const InaVecAVX& inVec2) {
        const __m256d testResult = _mm256_cmp_pd(inVec1.vec, inVec2.vec, _CMP_LE_OQ);
        const __m256d ones       = _mm256_set1_pd(1);
        return _mm256_and_pd(testResult, ones);
    }

    inline static InaVecAVX IsLower(const InaVecAVX& inVec1, const InaVecAVX& inVec2) {
        const __m256d testResult = _mm256_cmp_pd(inVec1.vec, inVec2.vec, _CMP_LT_OQ);
        const __m256d ones       = _mm256_set1_pd(1);
        return _mm256_and_pd(testResult, ones);
    }

    inline static InaVecAVX IsGreaterOrEqual(const InaVecAVX& inVec1, const InaVecAVX& inVec2) {
        const __m256d testResult = _mm256_cmp_pd(inVec1.vec, inVec2.vec, _CMP_GE_OQ);
        const __m256d ones       = _mm256_set1_pd(1);
        return _mm256_and_pd(testResult, ones);
    }

    inline static InaVecAVX IsGreater(const InaVecAVX& inVec1, const InaVecAVX& inVec2) {
        const __m256d testResult = _mm256_cmp_pd(inVec1.vec, inVec2.vec, _CMP_GT_OQ);
        const __m256d ones       = _mm256_set1_pd(1);
        return _mm256_and_pd(testResult, ones);
    }

    inline static InaVecAVX IsEqual(const InaVecAVX& inVec1, const InaVecAVX& inVec2) {
        const __m256d testResult = _mm256_cmp_pd(inVec1.vec, inVec2.vec, _CMP_EQ_OQ);
        const __m256d ones       = _mm256_set1_pd(1);
        return _mm256_and_pd(testResult, ones);
    }

    inline static InaVecAVX IsNotEqual(const InaVecAVX& inVec1, const InaVecAVX& inVec2) {
        const __m256d testResult = _mm256_cmp_pd(inVec1.vec, inVec2.vec, _CMP_NEQ_OQ);
        const __m256d ones       = _mm256_set1_pd(1);
        return _mm256_and_pd(testResult, ones);
    }

    inline static InaVecMaskAVX< double > IsLowerOrEqualMask(const InaVecAVX& inVec1, const InaVecAVX& inVec2) {
        return _mm256_castpd_si256(_mm256_cmp_pd(inVec1.vec, inVec2.vec, _CMP_LE_OQ));
    }

    inline static InaVecMaskAVX< double > IsLowerMask(const InaVecAVX& inVec1, const InaVecAVX& inVec2) {
        return _mm256_castpd_si256(_mm256_cmp_pd(inVec1.vec, inVec2.vec, _CMP_LT_OQ));
    }

    inline static InaVecMaskAVX< double > IsGreaterOrEqualMask(const InaVecAVX& inVec1, const InaVecAVX& inVec2) {
        return _mm256_castpd_si256(_mm256_cmp_pd(inVec1.vec, inVec2.vec, _CMP_GE_OQ));
    }

    inline static InaVecMaskAVX< double > IsGreaterMask(const InaVecAVX& inVec1, const InaVecAVX& inVec2) {
        return _mm256_castpd_si256(_mm256_cmp_pd(inVec1.vec, inVec2.vec, _CMP_GT_OQ));
    }

    inline static InaVecMaskAVX< double > IsEqualMask(const InaVecAVX& inVec1, const InaVecAVX& inVec2) {
        return _mm256_castpd_si256(_mm256_cmp_pd(inVec1.vec, inVec2.vec, _CMP_EQ_OQ));
    }

    inline static InaVecMaskAVX< double > IsNotEqualMask(const InaVecAVX& inVec1, const InaVecAVX& inVec2) {
        return _mm256_castpd_si256(_mm256_cmp_pd(inVec1.vec, inVec2.vec, _CMP_NEQ_OQ));
    }

    inline static InaVecAVX BitsAnd(const InaVecAVX& inVec1, const InaVecAVX& inVec2) {
        return _mm256_and_pd(inVec1.vec, inVec2.vec);
    }

    inline static InaVecAVX BitsNotAnd(const InaVecAVX& inVec1, const InaVecAVX& inVec2) {
        return _mm256_andnot_pd(inVec1.vec, inVec2.vec);
    }

    inline static InaVecAVX BitsOr(const InaVecAVX& inVec1, const InaVecAVX& inVec2) {
        return _mm256_or_pd(inVec1.vec, inVec2.vec);
    }

    inline static InaVecAVX BitsXor(const InaVecAVX& inVec1, const InaVecAVX& inVec2) {
        return _mm256_xor_pd(inVec1.vec, inVec2.vec);
    }

    inline static const char* GetName() {
        return "InaVecAVX<double>";
    }

    inline static InaIfElse< InaVecAVX< double > >::ThenClass If(const InaVecMaskAVX< double >& inTest) {
        return InaIfElse< InaVecAVX< double > >::IfClass().If(inTest);
    }

    inline static InaVecAVX IfElse(const InaVecMaskAVX< double >& inMask, const InaVecAVX& inIfTrue, const InaVecAVX& inIfFalse) {
        return _mm256_or_pd(IfTrue(inMask, inIfTrue.vec).vec,
                            IfFalse(inMask, inIfFalse.vec).vec);
    }

    inline static InaVecAVX IfTrue(const InaVecMaskAVX< double >& inMask, const InaVecAVX& inIfTrue) {
        return _mm256_and_pd(_mm256_castsi256_pd(inMask.getMask()), inIfTrue.vec);
    }

    inline static InaVecAVX IfFalse(const InaVecMaskAVX< double >& inMask, const InaVecAVX& inIfFalse) {
        return _mm256_andnot_pd(_mm256_castsi256_pd(inMask.getMask()), inIfFalse.vec);
    }

    // Inner operators
    inline InaVecAVX< double >& operator+=(const InaVecAVX< double >& inVec) {
        vec = _mm256_add_pd(vec, inVec.vec);
        return *this;
    }

    inline InaVecAVX< double >& operator-=(const InaVecAVX< double >& inVec) {
        vec = _mm256_sub_pd(vec, inVec.vec);
        return *this;
    }

    inline InaVecAVX< double >& operator/=(const InaVecAVX< double >& inVec) {
        vec = _mm256_div_pd(vec, inVec.vec);
        return *this;
    }

    inline InaVecAVX< double >& operator*=(const InaVecAVX< double >& inVec) {
        vec = _mm256_mul_pd(vec, inVec.vec);
        return *this;
    }

    inline InaVecAVX< double > operator-() const {
        const __m256d minus0 = _mm256_castsi256_pd(_mm256_set1_epi64x(static_cast< long long >(0x8000000000000000L)));
        return _mm256_xor_pd(vec, minus0);
    }

    inline InaVecAVX< double > pow(std::size_t power) const {
        return InaUtils::FastPow< InaVecAVX< double > >(*this, power);
    }
};

// Bits operators
inline InaVecAVX< double > operator&(const InaVecAVX< double >& inVec1, const InaVecAVX< double >& inVec2) {
    return InaVecAVX< double >::BitsAnd(inVec1, inVec2);
}

inline InaVecAVX< double > operator|(const InaVecAVX< double >& inVec1, const InaVecAVX< double >& inVec2) {
    return InaVecAVX< double >::BitsOr(inVec1, inVec2);
}

inline InaVecAVX< double > operator^(const InaVecAVX< double >& inVec1, const InaVecAVX< double >& inVec2) {
    return InaVecAVX< double >::BitsXor(inVec1, inVec2);
}

// Dual operators
inline InaVecAVX< double > operator+(const InaVecAVX< double >& inVec1, const InaVecAVX< double >& inVec2) {
    return _mm256_add_pd(inVec1.getVec(), inVec2.getVec());
}

inline InaVecAVX< double > operator-(const InaVecAVX< double >& inVec1, const InaVecAVX< double >& inVec2) {
    return _mm256_sub_pd(inVec1.getVec(), inVec2.getVec());
}

inline InaVecAVX< double > operator/(const InaVecAVX< double >& inVec1, const InaVecAVX< double >& inVec2) {
    return _mm256_div_pd(inVec1.getVec(), inVec2.getVec());
}

inline InaVecAVX< double > operator*(const InaVecAVX< double >& inVec1, const InaVecAVX< double >& inVec2) {
    return _mm256_mul_pd(inVec1.getVec(), inVec2.getVec());
}

// Tests and comparions
inline InaVecMaskAVX< double > operator<(const InaVecAVX< double >& inVec1, const InaVecAVX< double >& inVec2) {
    return InaVecAVX< double >::IsLowerMask(inVec1, inVec2);
}

inline InaVecMaskAVX< double > operator<=(const InaVecAVX< double >& inVec1, const InaVecAVX< double >& inVec2) {
    return InaVecAVX< double >::IsLowerOrEqualMask(inVec1, inVec2);
}

inline InaVecMaskAVX< double > operator>(const InaVecAVX< double >& inVec1, const InaVecAVX< double >& inVec2) {
    return InaVecAVX< double >::IsGreaterMask(inVec1, inVec2);
}

inline InaVecMaskAVX< double > operator>=(const InaVecAVX< double >& inVec1, const InaVecAVX< double >& inVec2) {
    return InaVecAVX< double >::IsGreaterOrEqualMask(inVec1, inVec2);
}

inline InaVecMaskAVX< double > operator==(const InaVecAVX< double >& inVec1, const InaVecAVX< double >& inVec2) {
    return InaVecAVX< double >::IsEqualMask(inVec1, inVec2);
}

inline InaVecMaskAVX< double > operator!=(const InaVecAVX< double >& inVec1, const InaVecAVX< double >& inVec2) {
    return InaVecAVX< double >::IsNotEqualMask(inVec1, inVec2);
}


#endif
