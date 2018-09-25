///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef INAVECAVX512COMMONDOUBLE_HPP
#define INAVECAVX512COMMONDOUBLE_HPP

#include "InastempConfig.h"
#include "Common/InaIfElse.hpp"
#include "Common/InaUtils.hpp"

#ifndef INASTEMP_USE_AVX512COMMON
#error InaVecAVX512COMMON<double> is included but AVX512COMMON is not enable in the configuration
#endif

#include "Common/InaFastExp.hpp"

#include <immintrin.h>

#include <cmath>
#include <initializer_list>

// Forward declarations
template < class RealType >
class InaVecMaskAVX512COMMON;

template < class RealType >
class InaVecAVX512COMMON;

// Mask type
template <>
class alignas(64) InaVecMaskAVX512COMMON< double > {
    __m512i mask;

public:
    // Classic constructors
    inline InaVecMaskAVX512COMMON() {
    }

    inline InaVecMaskAVX512COMMON(const InaVecMaskAVX512COMMON&) = default;
    inline InaVecMaskAVX512COMMON& operator=(const InaVecMaskAVX512COMMON&) = default;

    // Native data type compatibility
    inline /*not explicit*/ InaVecMaskAVX512COMMON(const __m512i inMask)
    : mask(inMask) {
    }

    inline InaVecMaskAVX512COMMON& operator=(const __m512i inMask) {
        mask = inMask;
        return (*this);
    }

    inline explicit operator __m512i() const {
        return mask;
    }

    inline __m512i getMask() const {
        return mask;
    }

    // Bool data type compatibility
    inline explicit InaVecMaskAVX512COMMON(const bool inBool) {
        mask = (inBool ? _mm512_set1_epi64(static_cast< long long >(0xFFFFFFFFFFFFFFFFL)) : _mm512_setzero_si512());
    }

    inline InaVecMaskAVX512COMMON& operator=(const bool inBool) {
        mask = (inBool ? _mm512_set1_epi64(static_cast< long long >(0xFFFFFFFFFFFFFFFFL)) : _mm512_setzero_si512());
        return (*this);
    }

    // Binary methods
    inline InaVecMaskAVX512COMMON Not() const {
        return NotAnd(mask, _mm512_set1_epi64(static_cast< long long >(0xFFFFFFFFFFFFFFFFL)));
    }

    inline bool isAllTrue() const {
        // true if all FF => !FF => 0 & FF => 0
        const __mmask8 testResult = _mm512_cmp_epu64_mask(mask, _mm512_set1_epi64(static_cast< long long >(0xFFFFFFFFFFFFFFFFL)), _MM_CMPINT_EQ);
        return testResult == 0xFF;
    }

    inline bool isAllFalse() const {
        // true if all zero
        const __mmask8 testResult = _mm512_cmp_epu64_mask(mask, _mm512_setzero_si512(), _MM_CMPINT_EQ);
        return testResult == 0xFF;
    }

    // Double args methods
    inline static InaVecMaskAVX512COMMON And(const InaVecMaskAVX512COMMON& inMask1, const InaVecMaskAVX512COMMON& inMask2) {
        return InaVecMaskAVX512COMMON(_mm512_and_si512(inMask1.mask, inMask2.mask));
    }

    inline static InaVecMaskAVX512COMMON NotAnd(const InaVecMaskAVX512COMMON& inMask1, const InaVecMaskAVX512COMMON& inMask2) {
        return InaVecMaskAVX512COMMON(_mm512_andnot_si512(inMask1.mask, inMask2.mask));
    }

    inline static InaVecMaskAVX512COMMON Or(const InaVecMaskAVX512COMMON& inMask1, const InaVecMaskAVX512COMMON& inMask2) {
        return InaVecMaskAVX512COMMON(_mm512_or_si512(inMask1.mask, inMask2.mask));
    }

    inline static InaVecMaskAVX512COMMON Xor(const InaVecMaskAVX512COMMON& inMask1, const InaVecMaskAVX512COMMON& inMask2) {
        return InaVecMaskAVX512COMMON(_mm512_xor_si512(inMask1.mask, inMask2.mask));
    }

    inline static bool IsEqual(const InaVecMaskAVX512COMMON& inMask1, const InaVecMaskAVX512COMMON& inMask2) {
        const __mmask8 testResult = _mm512_cmp_epu64_mask(inMask1.mask, inMask2.mask, _MM_CMPINT_EQ);
        return testResult == 0xFF;
    }

    inline static bool IsNotEqual(const InaVecMaskAVX512COMMON& inMask1, const InaVecMaskAVX512COMMON& inMask2) {
        const __mmask8 testResult = _mm512_cmp_epu64_mask(inMask1.mask, inMask2.mask, _MM_CMPINT_EQ);
        return testResult != 0xFF;
    }
};

// Mask must have operators
inline InaVecMaskAVX512COMMON< double > operator&(const InaVecMaskAVX512COMMON< double >& inMask1, const InaVecMaskAVX512COMMON< double >& inMask2) {
    return InaVecMaskAVX512COMMON< double >::And(inMask1, inMask2);
}

inline InaVecMaskAVX512COMMON< double > operator|(const InaVecMaskAVX512COMMON< double >& inMask1, const InaVecMaskAVX512COMMON< double >& inMask2) {
    return InaVecMaskAVX512COMMON< double >::Or(inMask1, inMask2);
}

inline InaVecMaskAVX512COMMON< double > operator^(const InaVecMaskAVX512COMMON< double >& inMask1, const InaVecMaskAVX512COMMON< double >& inMask2) {
    return InaVecMaskAVX512COMMON< double >::Xor(inMask1, inMask2);
}

inline bool operator==(const InaVecMaskAVX512COMMON< double >& inMask1, const InaVecMaskAVX512COMMON< double >& inMask2) {
    return InaVecMaskAVX512COMMON< double >::IsEqual(inMask1, inMask2);
}

inline bool operator!=(const InaVecMaskAVX512COMMON< double >& inMask1, const InaVecMaskAVX512COMMON< double >& inMask2) {
    return InaVecMaskAVX512COMMON< double >::IsNotEqual(inMask1, inMask2);
}

// Vec type
template <>
class alignas(64) InaVecAVX512COMMON< double > {
protected:
    __m512d vec;

public:
    using VecRawType            = __m512d;
    using MaskType              = InaVecMaskAVX512COMMON< double >;
    using RealType              = double;
    static const int VecLength  = 8;
    static const int Alignement = 64;

    inline InaVecAVX512COMMON() {
    }
    inline InaVecAVX512COMMON(const InaVecAVX512COMMON&) = default;
    inline InaVecAVX512COMMON& operator=(const InaVecAVX512COMMON&) = default;

    // Constructor from raw type
    inline /*not explicit*/ InaVecAVX512COMMON(const __m512d inVec)
    : vec(inVec) {
    }

    inline InaVecAVX512COMMON& operator=(const __m512d inVec) {
        vec = inVec;
        return *this;
    }

    inline void setFromRawType(const __m512d inVec) {
        vec = inVec;
    }

    inline explicit operator __m512d() const {
        return vec;
    }

    inline __m512d getVec() const {
        return vec;
    }

    // Constructor from scalar
    inline /*not explicit*/ InaVecAVX512COMMON(const double val)
    : vec(_mm512_set1_pd(val)) {
    }

    inline InaVecAVX512COMMON& operator=(const double val) {
        vec = _mm512_set1_pd(val);
        return *this;
    }

    inline void setFromScalar(const double val) {
        vec = _mm512_set1_pd(val);
    }

    // Constructor from vec
    inline InaVecAVX512COMMON(const std::initializer_list< double > lst)
    : InaVecAVX512COMMON(lst.begin()) {
    }

    inline explicit InaVecAVX512COMMON(const double ptr[])
    : vec(_mm512_loadu_pd(ptr)) {
    }

    inline InaVecAVX512COMMON& setFromArray(const double ptr[]) {
        vec = _mm512_loadu_pd(ptr);
        return *this;
    }

    inline InaVecAVX512COMMON& setFromAlignedArray(const double ptr[]) {
        vec = _mm512_load_pd(ptr);
        return *this;
    }

    inline InaVecAVX512COMMON& setFromIndirectArray(const double values[], const int inIndirection[]) {
        vec = _mm512_set_pd(
            values[inIndirection[7]],
            values[inIndirection[6]],
            values[inIndirection[5]],
            values[inIndirection[4]],
            values[inIndirection[3]],
            values[inIndirection[2]],
            values[inIndirection[1]],
            values[inIndirection[0]]);
        return *this;
    }

    inline InaVecAVX512COMMON& setFromIndirect2DArray(const double inArray[], const int inIndirection1[],
                                                      const int inLeadingDimension, const int inIndirection2[]) {
        vec = _mm512_set_pd(
            inArray[inIndirection1[7] * inLeadingDimension + inIndirection2[7]],
            inArray[inIndirection1[6] * inLeadingDimension + inIndirection2[6]],
            inArray[inIndirection1[5] * inLeadingDimension + inIndirection2[5]],
            inArray[inIndirection1[4] * inLeadingDimension + inIndirection2[4]],
            inArray[inIndirection1[3] * inLeadingDimension + inIndirection2[3]],
            inArray[inIndirection1[2] * inLeadingDimension + inIndirection2[2]],
            inArray[inIndirection1[1] * inLeadingDimension + inIndirection2[1]],
            inArray[inIndirection1[0] * inLeadingDimension + inIndirection2[0]]);
        return *this;
    }

    // Move back to array
    inline void storeInArray(double ptr[]) const {
        _mm512_storeu_pd(ptr, vec);
    }

    inline void storeInAlignedArray(double ptr[]) const {
        _mm512_store_pd(ptr, vec);
    }

    // Acce to individual values
    inline double at(const int index) const {
        alignas(Alignement) double allval[VecLength];
        _mm512_store_pd(allval, vec);
        return allval[index];
    }

    // Horizontal operation
    inline double horizontalSum() const {
#ifdef __INTEL_COMPILER
        return _mm512_reduce_add_pd(vec);
#else
        __m256d low  = _mm512_castpd512_pd256(vec);
        __m256d high = _mm512_extractf64x4_pd(vec, 1);
        __m256d val  = low + high;

        const __m128d valupper = _mm256_extractf128_pd(val, 1);
        const __m128d rest     = _mm256_castpd256_pd128(val);
        // Not in 512 _mm256_zeroupper();
        const __m128d valval = _mm_add_pd(valupper, rest);
        const __m128d res    = _mm_add_pd(_mm_permute_pd(valval, 1), valval);
        return _mm_cvtsd_f64(res);
#endif
    }

    inline double horizontalMul() const {
#ifdef __INTEL_COMPILER
        return _mm512_reduce_mul_pd(vec);
#else
        __m256d low  = _mm512_castpd512_pd256(vec);
        __m256d high = _mm512_extractf64x4_pd(vec, 1);
        __m256d val  = low * high;

        const __m128d valupper = _mm256_extractf128_pd(val, 1);
        const __m128d rest     = _mm256_castpd256_pd128(val);
        // Not in 512 _mm256_zeroupper();
        const __m128d valval = _mm_mul_pd(valupper, rest);
        const __m128d res    = _mm_mul_pd(_mm_permute_pd(valval, 1), valval);
        return _mm_cvtsd_f64(res);
#endif
    }

    inline InaVecAVX512COMMON sqrt() const {
        return _mm512_sqrt_pd(vec);
    }

    inline InaVecAVX512COMMON exp() const {
        const __m512d COEFF_LOG2E = _mm512_set1_pd(double(InaFastExp::CoeffLog2E()));
        const __m512d COEFF_A     = _mm512_set1_pd(double(InaFastExp::CoeffA64()));
        const __m512d COEFF_B     = _mm512_set1_pd(double(InaFastExp::CoeffB64()));
        const __m512d COEFF_P5_X  = _mm512_set1_pd(double(InaFastExp::GetCoefficient9_8()));
        const __m512d COEFF_P5_Y  = _mm512_set1_pd(double(InaFastExp::GetCoefficient9_7()));
        const __m512d COEFF_P5_Z  = _mm512_set1_pd(double(InaFastExp::GetCoefficient9_6()));
        const __m512d COEFF_P5_A  = _mm512_set1_pd(double(InaFastExp::GetCoefficient9_5()));
        const __m512d COEFF_P5_B  = _mm512_set1_pd(double(InaFastExp::GetCoefficient9_4()));
        const __m512d COEFF_P5_C  = _mm512_set1_pd(double(InaFastExp::GetCoefficient9_3()));
        const __m512d COEFF_P5_D  = _mm512_set1_pd(double(InaFastExp::GetCoefficient9_2()));
        const __m512d COEFF_P5_E  = _mm512_set1_pd(double(InaFastExp::GetCoefficient9_1()));
        const __m512d COEFF_P5_F  = _mm512_set1_pd(double(InaFastExp::GetCoefficient9_0()));

        __m512d x = _mm512_mul_pd(vec, COEFF_LOG2E);

        const __m512d fractional_part = _mm512_sub_pd(x, InaVecAVX512COMMON(x).floor().vec);

        __m512d factor = _mm512_add_pd(_mm512_mul_pd(_mm512_add_pd(
                                                         _mm512_mul_pd(_mm512_add_pd(_mm512_mul_pd(_mm512_add_pd(
                                                                                                       _mm512_mul_pd(_mm512_add_pd(_mm512_mul_pd(_mm512_add_pd(
                                                                                                                                                     _mm512_mul_pd(_mm512_add_pd(_mm512_mul_pd(_mm512_add_pd(_mm512_mul_pd(
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

        x = _mm512_sub_pd(x, factor);

        x = _mm512_add_pd(_mm512_mul_pd(COEFF_A, x), COEFF_B);

        alignas(64) double allvalreal[VecLength];
        _mm512_store_pd(allvalreal, x);

        alignas(64) long long int allvalint[VecLength] = { static_cast< long long int >(allvalreal[0]), static_cast< long long int >(allvalreal[1]),
                                                           static_cast< long long int >(allvalreal[2]), static_cast< long long int >(allvalreal[3]),
                                                           static_cast< long long int >(allvalreal[4]), static_cast< long long int >(allvalreal[5]),
                                                           static_cast< long long int >(allvalreal[6]), static_cast< long long int >(allvalreal[7]) };

        return _mm512_castsi512_pd(_mm512_load_epi64(reinterpret_cast< const __m512i* >(allvalint)));
    }

    inline InaVecAVX512COMMON expLowAcc() const {
        const __m512d COEFF_LOG2E = _mm512_set1_pd(double(InaFastExp::CoeffLog2E()));
        const __m512d COEFF_A     = _mm512_set1_pd(double(InaFastExp::CoeffA64()));
        const __m512d COEFF_B     = _mm512_set1_pd(double(InaFastExp::CoeffB64()));
        const __m512d COEFF_P5_C  = _mm512_set1_pd(double(InaFastExp::GetCoefficient4_3()));
        const __m512d COEFF_P5_D  = _mm512_set1_pd(double(InaFastExp::GetCoefficient4_2()));
        const __m512d COEFF_P5_E  = _mm512_set1_pd(double(InaFastExp::GetCoefficient4_1()));
        const __m512d COEFF_P5_F  = _mm512_set1_pd(double(InaFastExp::GetCoefficient4_0()));

        __m512d x = _mm512_mul_pd(vec, COEFF_LOG2E);

        const __m512d fractional_part = _mm512_sub_pd(x, InaVecAVX512COMMON(x).floor().vec);

        __m512d factor = _mm512_add_pd(_mm512_mul_pd(_mm512_add_pd(
                                                         _mm512_mul_pd(_mm512_add_pd(_mm512_mul_pd(
                                                                                         COEFF_P5_C, fractional_part),
                                                                                     COEFF_P5_D),
                                                                       fractional_part),
                                                         COEFF_P5_E),
                                                     fractional_part),
                                       COEFF_P5_F);

        x = _mm512_sub_pd(x, factor);

        x = _mm512_add_pd(_mm512_mul_pd(COEFF_A, x), COEFF_B);

        alignas(64) double allvalreal[VecLength];
        _mm512_store_pd(allvalreal, x);

        alignas(64) long long int allvalint[VecLength] = { static_cast< long long int >(allvalreal[0]), static_cast< long long int >(allvalreal[1]),
                                                           static_cast< long long int >(allvalreal[2]), static_cast< long long int >(allvalreal[3]),
                                                           static_cast< long long int >(allvalreal[4]), static_cast< long long int >(allvalreal[5]),
                                                           static_cast< long long int >(allvalreal[6]), static_cast< long long int >(allvalreal[7]) };

        return _mm512_castsi512_pd(_mm512_load_epi64(reinterpret_cast< const __m512i* >(allvalint)));
    }

    inline InaVecAVX512COMMON rsqrt() const {
        // _mm512_rsqrt28_pd(vec) => 1E-10 error
        return _mm512_div_pd(_mm512_set1_pd(1), _mm512_sqrt_pd(vec));
    }

    inline InaVecAVX512COMMON abs() const {
        const __m512d minus0 = _mm512_castsi512_pd(_mm512_set1_epi64(static_cast< long long >(0x8000000000000000L)));
        return _mm512_castsi512_pd(_mm512_andnot_epi32(_mm512_castpd_si512(minus0), _mm512_castpd_si512(vec)));
    }

    inline InaVecAVX512COMMON floor() const {
        return _mm512_cvt_roundps_pd(
            _mm256_cvtepi32_ps(
                _mm512_cvt_roundpd_epi32(vec, (_MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC))),
            _MM_FROUND_NO_EXC);
    }

    inline InaVecAVX512COMMON signOf() const {
        const __m512d minus0 = _mm512_castsi512_pd(_mm512_set1_epi64(static_cast< long long >(0x8000000000000000L)));
        const __m512d signs  = _mm512_castsi512_pd(_mm512_and_epi32(_mm512_castpd_si512(vec), _mm512_castpd_si512(minus0)));
        return _mm512_maskz_mov_pd(
            _mm512_cmp_pd_mask(_mm512_setzero_pd(), vec, _CMP_NEQ_OQ),
            _mm512_castsi512_pd(_mm512_or_epi32(_mm512_castpd_si512(signs),
                                                _mm512_castpd_si512(_mm512_set1_pd(1)))));
    }

    inline InaVecAVX512COMMON isPositive() const {
        const __mmask8 greater = _mm512_cmp_pd_mask(_mm512_setzero_pd(), vec, _CMP_LE_OQ);
        const __m512d ones     = _mm512_set1_pd(1);
        return _mm512_maskz_mov_pd(greater, ones);
    }

    inline InaVecAVX512COMMON isNegative() const {
        const __mmask8 less = _mm512_cmp_pd_mask(_mm512_setzero_pd(), vec, _CMP_GE_OQ);
        const __m512d ones  = _mm512_set1_pd(1);
        return _mm512_maskz_mov_pd(less, ones);
    }

    inline InaVecAVX512COMMON isPositiveStrict() const {
        const __mmask8 greater = _mm512_cmp_pd_mask(_mm512_setzero_pd(), vec, _CMP_LT_OQ);
        const __m512d ones     = _mm512_set1_pd(1);
        return _mm512_maskz_mov_pd(greater, ones);
    }

    inline InaVecAVX512COMMON isNegativeStrict() const {
        const __mmask8 less = _mm512_cmp_pd_mask(_mm512_setzero_pd(), vec, _CMP_GT_OQ);
        const __m512d ones  = _mm512_set1_pd(1);
        return _mm512_maskz_mov_pd(less, ones);
    }

    inline InaVecAVX512COMMON isZero() const {
        const __mmask8 equalZero = _mm512_cmp_pd_mask(_mm512_setzero_pd(), vec, _CMP_EQ_OQ);
        const __m512d ones       = _mm512_set1_pd(1);
        return _mm512_maskz_mov_pd(equalZero, ones);
    }

    inline InaVecAVX512COMMON isNotZero() const {
        const __mmask8 equalZero = _mm512_cmp_pd_mask(_mm512_setzero_pd(), vec, _CMP_NEQ_OQ);
        const __m512d ones       = _mm512_set1_pd(1);
        return _mm512_maskz_mov_pd(equalZero, ones);
    }

    inline InaVecMaskAVX512COMMON< double > isPositiveMask() const {
        return _mm512_castpd_si512(_mm512_maskz_mov_pd(_mm512_cmp_pd_mask(_mm512_setzero_pd(), vec, _CMP_LE_OQ),
                                                       _mm512_castsi512_pd(_mm512_set1_epi64(static_cast< long long >(0xFFFFFFFFFFFFFFFFL)))));
    }

    inline InaVecMaskAVX512COMMON< double > isNegativeMask() const {
        return _mm512_castpd_si512(_mm512_maskz_mov_pd(_mm512_cmp_pd_mask(_mm512_setzero_pd(), vec, _CMP_GE_OQ),
                                                       _mm512_castsi512_pd(_mm512_set1_epi64(static_cast< long long >(0xFFFFFFFFFFFFFFFFL)))));
    }

    inline InaVecMaskAVX512COMMON< double > isPositiveStrictMask() const {
        return _mm512_castpd_si512(_mm512_maskz_mov_pd(_mm512_cmp_pd_mask(_mm512_setzero_pd(), vec, _CMP_LT_OQ),
                                                       _mm512_castsi512_pd(_mm512_set1_epi64(static_cast< long long >(0xFFFFFFFFFFFFFFFFL)))));
    }

    inline InaVecMaskAVX512COMMON< double > isNegativeStrictMask() const {
        return _mm512_castpd_si512(_mm512_maskz_mov_pd(_mm512_cmp_pd_mask(_mm512_setzero_pd(), vec, _CMP_GT_OQ),
                                                       _mm512_castsi512_pd(_mm512_set1_epi64(static_cast< long long >(0xFFFFFFFFFFFFFFFFL)))));
    }

    inline InaVecMaskAVX512COMMON< double > isZeroMask() const {
        return _mm512_castpd_si512(_mm512_maskz_mov_pd(_mm512_cmp_pd_mask(_mm512_setzero_pd(), vec, _CMP_EQ_OQ),
                                                       _mm512_castsi512_pd(_mm512_set1_epi64(static_cast< long long >(0xFFFFFFFFFFFFFFFFL)))));
    }

    inline InaVecMaskAVX512COMMON< double > isNotZeroMask() const {
        return _mm512_castpd_si512(_mm512_maskz_mov_pd(_mm512_cmp_pd_mask(_mm512_setzero_pd(), vec, _CMP_NEQ_OQ),
                                                       _mm512_castsi512_pd(_mm512_set1_epi64(static_cast< long long >(0xFFFFFFFFFFFFFFFFL)))));
    }

    // Static basic methods
    inline static InaVecAVX512COMMON GetZero() {
        return InaVecAVX512COMMON(_mm512_setzero_pd());
    }

    inline static InaVecAVX512COMMON GetOne() {
        return InaVecAVX512COMMON(_mm512_set1_pd(1));
    }

    inline static InaVecAVX512COMMON Min(const InaVecAVX512COMMON& inVec1, const InaVecAVX512COMMON& inVec2) {
        return _mm512_min_pd(inVec1.vec, inVec2.vec);
    }

    inline static InaVecAVX512COMMON Max(const InaVecAVX512COMMON& inVec1, const InaVecAVX512COMMON& inVec2) {
        return _mm512_max_pd(inVec1.vec, inVec2.vec);
    }

    inline static InaVecAVX512COMMON IsLowerOrEqual(const InaVecAVX512COMMON& inVec1, const InaVecAVX512COMMON& inVec2) {
        const __mmask8 testResult = _mm512_cmp_pd_mask(inVec1.vec, inVec2.vec, _CMP_LE_OQ);
        const __m512d ones        = _mm512_set1_pd(1);
        return _mm512_maskz_mov_pd(testResult, ones);
    }

    inline static InaVecAVX512COMMON IsLower(const InaVecAVX512COMMON& inVec1, const InaVecAVX512COMMON& inVec2) {
        const __mmask8 testResult = _mm512_cmp_pd_mask(inVec1.vec, inVec2.vec, _CMP_LT_OQ);
        const __m512d ones        = _mm512_set1_pd(1);
        return _mm512_maskz_mov_pd(testResult, ones);
    }

    inline static InaVecAVX512COMMON IsGreaterOrEqual(const InaVecAVX512COMMON& inVec1, const InaVecAVX512COMMON& inVec2) {
        const __mmask8 testResult = _mm512_cmp_pd_mask(inVec1.vec, inVec2.vec, _CMP_GE_OQ);
        const __m512d ones        = _mm512_set1_pd(1);
        return _mm512_maskz_mov_pd(testResult, ones);
    }

    inline static InaVecAVX512COMMON IsGreater(const InaVecAVX512COMMON& inVec1, const InaVecAVX512COMMON& inVec2) {
        const __mmask8 testResult = _mm512_cmp_pd_mask(inVec1.vec, inVec2.vec, _CMP_GT_OQ);
        const __m512d ones        = _mm512_set1_pd(1);
        return _mm512_maskz_mov_pd(testResult, ones);
    }

    inline static InaVecAVX512COMMON IsEqual(const InaVecAVX512COMMON& inVec1, const InaVecAVX512COMMON& inVec2) {
        const __mmask8 testResult = _mm512_cmp_pd_mask(inVec1.vec, inVec2.vec, _CMP_EQ_OQ);
        const __m512d ones        = _mm512_set1_pd(1);
        return _mm512_maskz_mov_pd(testResult, ones);
    }

    inline static InaVecAVX512COMMON IsNotEqual(const InaVecAVX512COMMON& inVec1, const InaVecAVX512COMMON& inVec2) {
        const __mmask8 testResult = _mm512_cmp_pd_mask(inVec1.vec, inVec2.vec, _CMP_NEQ_OQ);
        const __m512d ones        = _mm512_set1_pd(1);
        return _mm512_maskz_mov_pd(testResult, ones);
    }

    inline static InaVecMaskAVX512COMMON< double > IsLowerOrEqualMask(const InaVecAVX512COMMON& inVec1, const InaVecAVX512COMMON& inVec2) {
        return _mm512_castpd_si512(_mm512_maskz_mov_pd(_mm512_cmp_pd_mask(inVec1.vec, inVec2.vec, _CMP_LE_OQ),
                                                       _mm512_castsi512_pd(_mm512_set1_epi64(static_cast< long long >(0xFFFFFFFFFFFFFFFFL)))));
    }

    inline static InaVecMaskAVX512COMMON< double > IsLowerMask(const InaVecAVX512COMMON& inVec1, const InaVecAVX512COMMON& inVec2) {
        return _mm512_castpd_si512(_mm512_maskz_mov_pd(_mm512_cmp_pd_mask(inVec1.vec, inVec2.vec, _CMP_LT_OQ),
                                                       _mm512_castsi512_pd(_mm512_set1_epi64(static_cast< long long >(0xFFFFFFFFFFFFFFFFL)))));
    }

    inline static InaVecMaskAVX512COMMON< double > IsGreaterOrEqualMask(const InaVecAVX512COMMON& inVec1, const InaVecAVX512COMMON& inVec2) {
        return _mm512_castpd_si512(_mm512_maskz_mov_pd(_mm512_cmp_pd_mask(inVec1.vec, inVec2.vec, _CMP_GE_OQ),
                                                       _mm512_castsi512_pd(_mm512_set1_epi64(static_cast< long long >(0xFFFFFFFFFFFFFFFFL)))));
    }

    inline static InaVecMaskAVX512COMMON< double > IsGreaterMask(const InaVecAVX512COMMON& inVec1, const InaVecAVX512COMMON& inVec2) {
        return _mm512_castpd_si512(_mm512_maskz_mov_pd(_mm512_cmp_pd_mask(inVec1.vec, inVec2.vec, _CMP_GT_OQ),
                                                       _mm512_castsi512_pd(_mm512_set1_epi64(static_cast< long long >(0xFFFFFFFFFFFFFFFFL)))));
    }

    inline static InaVecMaskAVX512COMMON< double > IsEqualMask(const InaVecAVX512COMMON& inVec1, const InaVecAVX512COMMON& inVec2) {
        return _mm512_castpd_si512(_mm512_maskz_mov_pd(_mm512_cmp_pd_mask(inVec1.vec, inVec2.vec, _CMP_EQ_OQ),
                                                       _mm512_castsi512_pd(_mm512_set1_epi64(static_cast< long long >(0xFFFFFFFFFFFFFFFFL)))));
    }

    inline static InaVecMaskAVX512COMMON< double > IsNotEqualMask(const InaVecAVX512COMMON& inVec1, const InaVecAVX512COMMON& inVec2) {
        return _mm512_castpd_si512(_mm512_maskz_mov_pd(_mm512_cmp_pd_mask(inVec1.vec, inVec2.vec, _CMP_NEQ_OQ),
                                                       _mm512_castsi512_pd(_mm512_set1_epi64(static_cast< long long >(0xFFFFFFFFFFFFFFFFL)))));
    }

    inline static InaVecAVX512COMMON BitsAnd(const InaVecAVX512COMMON& inVec1, const InaVecAVX512COMMON& inVec2) {
        return _mm512_castsi512_pd(_mm512_and_si512(_mm512_castpd_si512(inVec1.vec), _mm512_castpd_si512(inVec2.vec)));
    }

    inline static InaVecAVX512COMMON BitsNotAnd(const InaVecAVX512COMMON& inVec1, const InaVecAVX512COMMON& inVec2) {
        return _mm512_castsi512_pd(_mm512_andnot_si512(_mm512_castpd_si512(inVec1.vec), _mm512_castpd_si512(inVec2.vec)));
    }

    inline static InaVecAVX512COMMON BitsOr(const InaVecAVX512COMMON& inVec1, const InaVecAVX512COMMON& inVec2) {
        return _mm512_castsi512_pd(_mm512_or_si512(_mm512_castpd_si512(inVec1.vec), _mm512_castpd_si512(inVec2.vec)));
    }

    inline static InaVecAVX512COMMON BitsXor(const InaVecAVX512COMMON& inVec1, const InaVecAVX512COMMON& inVec2) {
        return _mm512_castsi512_pd(_mm512_xor_si512(_mm512_castpd_si512(inVec1.vec), _mm512_castpd_si512(inVec2.vec)));
    }

    inline static const char* GetName() {
        return "InaVecAVX512COMMON<double>";
    }

    inline static InaIfElse< InaVecAVX512COMMON< double > >::ThenClass If(const InaVecMaskAVX512COMMON< double >& inTest) {
        return InaIfElse< InaVecAVX512COMMON< double > >::IfClass().If(inTest);
    }

    inline static InaVecAVX512COMMON IfElse(const InaVecMaskAVX512COMMON< double >& inMask, const InaVecAVX512COMMON& inIfTrue, const InaVecAVX512COMMON& inIfFalse) {
        return _mm512_castsi512_pd(_mm512_or_si512(_mm512_castpd_si512(IfTrue(inMask, inIfTrue.vec).vec),
                                                   _mm512_castpd_si512(IfFalse(inMask, inIfFalse.vec).vec)));
    }

    inline static InaVecAVX512COMMON IfTrue(const InaVecMaskAVX512COMMON< double >& inMask, const InaVecAVX512COMMON& inIfTrue) {
        return _mm512_castsi512_pd(_mm512_and_si512(inMask.getMask(), _mm512_castpd_si512(inIfTrue.vec)));
    }

    inline static InaVecAVX512COMMON IfFalse(const InaVecMaskAVX512COMMON< double >& inMask, const InaVecAVX512COMMON& inIfFalse) {
        return _mm512_castsi512_pd(_mm512_andnot_si512(inMask.getMask(), _mm512_castpd_si512(inIfFalse.vec)));
    }

    // Inner operators
    inline InaVecAVX512COMMON< double >& operator+=(const InaVecAVX512COMMON< double >& inVec) {
        vec = _mm512_add_pd(vec, inVec.vec);
        return *this;
    }

    inline InaVecAVX512COMMON< double >& operator-=(const InaVecAVX512COMMON< double >& inVec) {
        vec = _mm512_sub_pd(vec, inVec.vec);
        return *this;
    }

    inline InaVecAVX512COMMON< double >& operator/=(const InaVecAVX512COMMON< double >& inVec) {
        vec = _mm512_div_pd(vec, inVec.vec);
        return *this;
    }

    inline InaVecAVX512COMMON< double >& operator*=(const InaVecAVX512COMMON< double >& inVec) {
        vec = _mm512_mul_pd(vec, inVec.vec);
        return *this;
    }

    inline InaVecAVX512COMMON< double > operator-() const {
        const __m512d minus0 = _mm512_castsi512_pd(_mm512_set1_epi64(static_cast< long long >(0x8000000000000000L)));
        return _mm512_castsi512_pd(_mm512_xor_si512(_mm512_castpd_si512(vec), _mm512_castpd_si512(minus0)));
    }

    inline InaVecAVX512COMMON< double > pow(std::size_t power) const {
        return InaUtils::FastPow< InaVecAVX512COMMON< double > >(*this, power);
    }
};

// Bits operators
inline InaVecAVX512COMMON< double > operator&(const InaVecAVX512COMMON< double >& inVec1, const InaVecAVX512COMMON< double >& inVec2) {
    return InaVecAVX512COMMON< double >::BitsAnd(inVec1, inVec2);
}

inline InaVecAVX512COMMON< double > operator|(const InaVecAVX512COMMON< double >& inVec1, const InaVecAVX512COMMON< double >& inVec2) {
    return InaVecAVX512COMMON< double >::BitsOr(inVec1, inVec2);
}

inline InaVecAVX512COMMON< double > operator^(const InaVecAVX512COMMON< double >& inVec1, const InaVecAVX512COMMON< double >& inVec2) {
    return InaVecAVX512COMMON< double >::BitsXor(inVec1, inVec2);
}

// Dual operators
inline InaVecAVX512COMMON< double > operator+(const InaVecAVX512COMMON< double >& inVec1, const InaVecAVX512COMMON< double >& inVec2) {
    return _mm512_add_pd(inVec1.getVec(), inVec2.getVec());
}

inline InaVecAVX512COMMON< double > operator-(const InaVecAVX512COMMON< double >& inVec1, const InaVecAVX512COMMON< double >& inVec2) {
    return _mm512_sub_pd(inVec1.getVec(), inVec2.getVec());
}

inline InaVecAVX512COMMON< double > operator/(const InaVecAVX512COMMON< double >& inVec1, const InaVecAVX512COMMON< double >& inVec2) {
    return _mm512_div_pd(inVec1.getVec(), inVec2.getVec());
}

inline InaVecAVX512COMMON< double > operator*(const InaVecAVX512COMMON< double >& inVec1, const InaVecAVX512COMMON< double >& inVec2) {
    return _mm512_mul_pd(inVec1.getVec(), inVec2.getVec());
}

// Tests and comparions
inline InaVecMaskAVX512COMMON< double > operator<(const InaVecAVX512COMMON< double >& inVec1, const InaVecAVX512COMMON< double >& inVec2) {
    return InaVecAVX512COMMON< double >::IsLowerMask(inVec1, inVec2);
}

inline InaVecMaskAVX512COMMON< double > operator<=(const InaVecAVX512COMMON< double >& inVec1, const InaVecAVX512COMMON< double >& inVec2) {
    return InaVecAVX512COMMON< double >::IsLowerOrEqualMask(inVec1, inVec2);
}

inline InaVecMaskAVX512COMMON< double > operator>(const InaVecAVX512COMMON< double >& inVec1, const InaVecAVX512COMMON< double >& inVec2) {
    return InaVecAVX512COMMON< double >::IsGreaterMask(inVec1, inVec2);
}

inline InaVecMaskAVX512COMMON< double > operator>=(const InaVecAVX512COMMON< double >& inVec1, const InaVecAVX512COMMON< double >& inVec2) {
    return InaVecAVX512COMMON< double >::IsGreaterOrEqualMask(inVec1, inVec2);
}

inline InaVecMaskAVX512COMMON< double > operator==(const InaVecAVX512COMMON< double >& inVec1, const InaVecAVX512COMMON< double >& inVec2) {
    return InaVecAVX512COMMON< double >::IsEqualMask(inVec1, inVec2);
}

inline InaVecMaskAVX512COMMON< double > operator!=(const InaVecAVX512COMMON< double >& inVec1, const InaVecAVX512COMMON< double >& inVec2) {
    return InaVecAVX512COMMON< double >::IsNotEqualMask(inVec1, inVec2);
}


#endif
