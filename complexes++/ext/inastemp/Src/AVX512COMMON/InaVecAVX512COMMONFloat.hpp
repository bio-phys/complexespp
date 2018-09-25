///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef INAVECAVX512COMMONFLOAT_HPP
#define INAVECAVX512COMMONFLOAT_HPP

#include "InastempConfig.h"
#include "Common/InaIfElse.hpp"
#include "Common/InaUtils.hpp"

#ifndef INASTEMP_USE_AVX512COMMON
#error InaVecAVX512COMMON512COMMON<float> is included but AVX512COMMON is not enable in the configuration
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
class alignas(64) InaVecMaskAVX512COMMON< float > {
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
        mask = (inBool ? _mm512_set1_epi32(static_cast< int >(0xFFFFFFFF)) : _mm512_setzero_si512());
    }

    inline InaVecMaskAVX512COMMON& operator=(const bool inBool) {
        mask = (inBool ? _mm512_set1_epi32(static_cast< int >(0xFFFFFFFF)) : _mm512_setzero_si512());
        return (*this);
    }

    // Binary methods
    inline InaVecMaskAVX512COMMON Not() const {
        return NotAnd(mask, _mm512_set1_epi32(static_cast< int >(0xFFFFFFFF)));
    }


    inline bool isAllTrue() const {
        // true if all FF => !FF => 0 & FF => 0
        // sde bug here with _mm512_cmp_epu32_mask return 0xFF instead of 0xFFFF so use 64 bits version
        const __mmask8 testResult = _mm512_cmp_epu64_mask(mask, _mm512_set1_epi64(0xFFFFFFFFFFFFFFFFUL), _MM_CMPINT_EQ);
        return testResult == 0xFF;
    }

    inline bool isAllFalse() const {
        // true if all zero
        // sde bug here with _mm512_cmp_epu32_mask return 0xFF instead of 0xFFFF so use 64 bits version
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
        // sde bug here with _mm512_cmp_epu32_mask return 0xFF instead of 0xFFFF so use 64 bits version
        const __mmask8 testResult = _mm512_cmp_epu64_mask(inMask1.mask, inMask2.mask, _MM_CMPINT_EQ);
        return testResult == 0xFF;
    }

    inline static bool IsNotEqual(const InaVecMaskAVX512COMMON& inMask1, const InaVecMaskAVX512COMMON& inMask2) {
        // sde bug here with _mm512_cmp_epu32_mask return 0xFF instead of 0xFFFF so use 64 bits version
        const __mmask8 testResult = _mm512_cmp_epu64_mask(inMask1.mask, inMask2.mask, _MM_CMPINT_EQ);
        return testResult != 0xFF;
    }
};

// Mask must have operators
inline InaVecMaskAVX512COMMON< float > operator&(const InaVecMaskAVX512COMMON< float >& inMask1, const InaVecMaskAVX512COMMON< float >& inMask2) {
    return InaVecMaskAVX512COMMON< float >::And(inMask1, inMask2);
}

inline InaVecMaskAVX512COMMON< float > operator|(const InaVecMaskAVX512COMMON< float >& inMask1, const InaVecMaskAVX512COMMON< float >& inMask2) {
    return InaVecMaskAVX512COMMON< float >::Or(inMask1, inMask2);
}

inline InaVecMaskAVX512COMMON< float > operator^(const InaVecMaskAVX512COMMON< float >& inMask1, const InaVecMaskAVX512COMMON< float >& inMask2) {
    return InaVecMaskAVX512COMMON< float >::Xor(inMask1, inMask2);
}

inline bool operator==(const InaVecMaskAVX512COMMON< float >& inMask1, const InaVecMaskAVX512COMMON< float >& inMask2) {
    return InaVecMaskAVX512COMMON< float >::IsEqual(inMask1, inMask2);
}

inline bool operator!=(const InaVecMaskAVX512COMMON< float >& inMask1, const InaVecMaskAVX512COMMON< float >& inMask2) {
    return InaVecMaskAVX512COMMON< float >::IsNotEqual(inMask1, inMask2);
}

// Vec type
template <>
class alignas(64) InaVecAVX512COMMON< float > {
protected:
    __m512 vec;

public:
    using VecRawType            = __m512;
    using MaskType              = InaVecMaskAVX512COMMON< float >;
    using RealType              = float;
    static const int VecLength  = 16;
    static const int Alignement = 64;

    inline InaVecAVX512COMMON() {
    }
    inline InaVecAVX512COMMON(const InaVecAVX512COMMON&) = default;
    inline InaVecAVX512COMMON& operator=(const InaVecAVX512COMMON&) = default;

    // Constructor from raw type
    inline /*not explicit*/ InaVecAVX512COMMON(const __m512 inVec)
    : vec(inVec) {
    }

    inline InaVecAVX512COMMON& operator=(const __m512 inVec) {
        vec = inVec;
        return *this;
    }

    inline void setFromRawType(const __m512 inVec) {
        vec = inVec;
    }

    inline explicit operator __m512() const {
        return vec;
    }

    inline __m512 getVec() const {
        return vec;
    }

    // Constructor from scalar
    inline /*not explicit*/ InaVecAVX512COMMON(const float val)
    : vec(_mm512_set1_ps(val)) {
    }

    inline InaVecAVX512COMMON& operator=(const float val) {
        vec = _mm512_set1_ps(val);
        return *this;
    }

    inline void setFromScalar(const float val) {
        vec = _mm512_set1_ps(val);
    }

    // Constructor from vec
    inline InaVecAVX512COMMON(const std::initializer_list< float > lst)
    : InaVecAVX512COMMON(lst.begin()) {
    }

    inline explicit InaVecAVX512COMMON(const float ptr[])
    : vec(_mm512_loadu_ps(ptr)) {
    }

    inline InaVecAVX512COMMON& setFromArray(const float ptr[]) {
        vec = _mm512_loadu_ps(ptr);
        return *this;
    }

    inline InaVecAVX512COMMON& setFromAlignedArray(const float ptr[]) {
        vec = _mm512_load_ps(ptr);
        return *this;
    }

    inline InaVecAVX512COMMON& setFromIndirectArray(const float values[], const int inIndirection[]) {
        vec = _mm512_set_ps(
            values[inIndirection[15]],
            values[inIndirection[14]],
            values[inIndirection[13]],
            values[inIndirection[12]],
            values[inIndirection[11]],
            values[inIndirection[10]],
            values[inIndirection[9]],
            values[inIndirection[8]],
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

    inline InaVecAVX512COMMON& setFromIndirect2DArray(const float inArray[], const int inIndirection1[],
                                                      const int inLeadingDimension, const int inIndirection2[]) {
        vec = _mm512_set_ps(
            inArray[inIndirection1[15] * inLeadingDimension + inIndirection2[15]],
            inArray[inIndirection1[14] * inLeadingDimension + inIndirection2[14]],
            inArray[inIndirection1[13] * inLeadingDimension + inIndirection2[13]],
            inArray[inIndirection1[12] * inLeadingDimension + inIndirection2[12]],
            inArray[inIndirection1[11] * inLeadingDimension + inIndirection2[11]],
            inArray[inIndirection1[10] * inLeadingDimension + inIndirection2[10]],
            inArray[inIndirection1[9] * inLeadingDimension + inIndirection2[9]],
            inArray[inIndirection1[8] * inLeadingDimension + inIndirection2[8]],
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
    inline void storeInArray(float ptr[]) const {
        _mm512_storeu_ps(ptr, vec);
    }

    inline void storeInAlignedArray(float ptr[]) const {
        _mm512_store_ps(ptr, vec);
    }

    // Acce to individual values
    inline float at(const int index) const {
        alignas(Alignement) float allval[VecLength];
        _mm512_store_ps(allval, vec);
        return allval[index];
    }

    // Horizontal operation
    inline float horizontalSum() const {
#ifdef __INTEL_COMPILER
        return _mm512_reduce_add_ps(vec);
#else
        __m256 low  = _mm512_castps512_ps256(vec);
        __m256 high = _mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(vec), 1));
        __m256 val  = low + high;

        const __m128 valupper = _mm256_extractf128_ps(val, 1);
        const __m128 rest     = _mm256_extractf128_ps(val, 0);
        // Not in 512 _mm256_zeroupper();
        const __m128 valval = _mm_add_ps(valupper,
                                         rest);
        __m128 valsum = _mm_add_ps(_mm_permute_ps(valval, 0x1B), valval);
        __m128 res    = _mm_add_ps(_mm_permute_ps(valsum, 0xB1), valsum);
        return _mm_cvtss_f32(res);
#endif
    }

    inline float horizontalMul() const {
#ifdef __INTEL_COMPILER
        return _mm512_reduce_mul_ps(vec);
#else
        __m256 low  = _mm512_castps512_ps256(vec);
        __m256 high = _mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(vec), 1));
        __m256 val  = low * high;

        const __m128 valupper = _mm256_extractf128_ps(val, 1);
        const __m128 rest     = _mm256_extractf128_ps(val, 0);
        // Not in 512 _mm256_zeroupper();
        const __m128 valval = _mm_mul_ps(valupper,
                                         rest);
        __m128 valsum = _mm_mul_ps(_mm_permute_ps(valval, 0x1B), valval);
        __m128 res    = _mm_mul_ps(_mm_permute_ps(valsum, 0xB1), valsum);
        return _mm_cvtss_f32(res);
#endif
    }

    inline InaVecAVX512COMMON sqrt() const {
        return _mm512_sqrt_ps(vec);
    }

    inline InaVecAVX512COMMON exp() const {
        const __m512 COEFF_LOG2E = _mm512_set1_ps(float(InaFastExp::CoeffLog2E()));
        const __m512 COEFF_A     = _mm512_set1_ps(float(InaFastExp::CoeffA32()));
        const __m512 COEFF_B     = _mm512_set1_ps(float(InaFastExp::CoeffB32()));
        const __m512 COEFF_P5_A  = _mm512_set1_ps(float(InaFastExp::GetCoefficient6_5()));
        const __m512 COEFF_P5_B  = _mm512_set1_ps(float(InaFastExp::GetCoefficient6_4()));
        const __m512 COEFF_P5_C  = _mm512_set1_ps(float(InaFastExp::GetCoefficient6_3()));
        const __m512 COEFF_P5_D  = _mm512_set1_ps(float(InaFastExp::GetCoefficient6_2()));
        const __m512 COEFF_P5_E  = _mm512_set1_ps(float(InaFastExp::GetCoefficient6_1()));
        const __m512 COEFF_P5_F  = _mm512_set1_ps(float(InaFastExp::GetCoefficient6_0()));

        __m512 x = _mm512_mul_ps(vec, COEFF_LOG2E);

        const __m512 fractional_part = _mm512_sub_ps(x, InaVecAVX512COMMON(x).floor().vec);

        __m512 factor = _mm512_add_ps(_mm512_mul_ps(_mm512_add_ps(_mm512_mul_ps(_mm512_add_ps(
                                                                                    _mm512_mul_ps(_mm512_add_ps(_mm512_mul_ps(_mm512_add_ps(_mm512_mul_ps(
                                                                                                                                                COEFF_P5_A, fractional_part),
                                                                                                                                            COEFF_P5_B),
                                                                                                                              fractional_part),
                                                                                                                COEFF_P5_C),
                                                                                                  fractional_part),
                                                                                    COEFF_P5_D),
                                                                                fractional_part),
                                                                  COEFF_P5_E),
                                                    fractional_part),
                                      COEFF_P5_F);

        x = _mm512_sub_ps(x, factor);

        __m512i castedInteger = _mm512_cvtps_epi32(_mm512_add_ps(_mm512_mul_ps(COEFF_A, x), COEFF_B));

        return _mm512_castsi512_ps(castedInteger);
    }

    inline InaVecAVX512COMMON expLowAcc() const {
        const __m512 COEFF_LOG2E = _mm512_set1_ps(float(InaFastExp::CoeffLog2E()));
        const __m512 COEFF_A     = _mm512_set1_ps(float(InaFastExp::CoeffA32()));
        const __m512 COEFF_B     = _mm512_set1_ps(float(InaFastExp::CoeffB32()));
        const __m512 COEFF_P5_D  = _mm512_set1_ps(float(InaFastExp::GetCoefficient3_2()));
        const __m512 COEFF_P5_E  = _mm512_set1_ps(float(InaFastExp::GetCoefficient3_1()));
        const __m512 COEFF_P5_F  = _mm512_set1_ps(float(InaFastExp::GetCoefficient3_0()));

        __m512 x = _mm512_mul_ps(vec, COEFF_LOG2E);

        const __m512 fractional_part = _mm512_sub_ps(x, InaVecAVX512COMMON(x).floor().vec);

        __m512 factor = _mm512_add_ps(_mm512_mul_ps(
                                          _mm512_add_ps(_mm512_mul_ps(
                                                            COEFF_P5_D, fractional_part),
                                                        COEFF_P5_E),
                                          fractional_part),
                                      COEFF_P5_F);

        x = _mm512_sub_ps(x, factor);

        __m512i castedInteger = _mm512_cvtps_epi32(_mm512_add_ps(_mm512_mul_ps(COEFF_A, x), COEFF_B));

        return _mm512_castsi512_ps(castedInteger);
    }

    inline InaVecAVX512COMMON rsqrt() const {
        return _mm512_rsqrt28_ps(vec);
    }

    inline InaVecAVX512COMMON abs() const {
        const __m512 minus0 = _mm512_castsi512_ps(_mm512_set1_epi32(static_cast< int >(0x80000000)));
        return _mm512_castsi512_ps(_mm512_andnot_epi32(_mm512_castps_si512(minus0), _mm512_castps_si512(vec)));
    }

    inline InaVecAVX512COMMON floor() const {
        return _mm512_cvt_roundepi32_ps(
            _mm512_cvt_roundps_epi32(vec, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)),
            (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
    }

    inline InaVecAVX512COMMON signOf() const {
        const __m512 minus0 = _mm512_castsi512_ps(_mm512_set1_epi32(static_cast< int >(0x80000000)));
        const __m512 signs  = _mm512_castsi512_ps(_mm512_and_epi32(_mm512_castps_si512(vec), _mm512_castps_si512(minus0)));
        return _mm512_maskz_mov_ps(
            _mm512_cmp_ps_mask(_mm512_setzero_ps(), vec, _CMP_NEQ_OQ),
            _mm512_castsi512_ps(_mm512_or_epi32(_mm512_castps_si512(signs),
                                                _mm512_castps_si512(_mm512_set1_ps(1)))));
    }

    inline InaVecAVX512COMMON isPositive() const {
        const __mmask16 greater = _mm512_cmp_ps_mask(_mm512_setzero_ps(), vec, _CMP_LE_OQ);
        const __m512 ones       = _mm512_set1_ps(1);
        return _mm512_maskz_mov_ps(greater, ones);
    }

    inline InaVecAVX512COMMON isNegative() const {
        const __mmask16 less = _mm512_cmp_ps_mask(_mm512_setzero_ps(), vec, _CMP_GE_OQ);
        const __m512 ones    = _mm512_set1_ps(1);
        return _mm512_maskz_mov_ps(less, ones);
    }

    inline InaVecAVX512COMMON isPositiveStrict() const {
        const __mmask16 greater = _mm512_cmp_ps_mask(_mm512_setzero_ps(), vec, _CMP_LT_OQ);
        const __m512 ones       = _mm512_set1_ps(1);
        return _mm512_maskz_mov_ps(greater, ones);
    }

    inline InaVecAVX512COMMON isNegativeStrict() const {
        const __mmask16 less = _mm512_cmp_ps_mask(_mm512_setzero_ps(), vec, _CMP_GT_OQ);
        const __m512 ones    = _mm512_set1_ps(1);
        return _mm512_maskz_mov_ps(less, ones);
    }

    inline InaVecAVX512COMMON isZero() const {
        const __mmask16 equalZero = _mm512_cmp_ps_mask(_mm512_setzero_ps(), vec, _CMP_EQ_OQ);
        const __m512 ones         = _mm512_set1_ps(1);
        return _mm512_maskz_mov_ps(equalZero, ones);
    }

    inline InaVecAVX512COMMON isNotZero() const {
        const __mmask16 equalZero = _mm512_cmp_ps_mask(_mm512_setzero_ps(), vec, _CMP_NEQ_OQ);
        const __m512 ones         = _mm512_set1_ps(1);
        return _mm512_maskz_mov_ps(equalZero, ones);
    }

    inline InaVecMaskAVX512COMMON< float > isPositiveMask() const {
        return _mm512_castps_si512(_mm512_maskz_mov_ps(_mm512_cmp_ps_mask(_mm512_setzero_ps(), vec, _CMP_LE_OQ),
                                                       _mm512_castsi512_ps(_mm512_set1_epi32(static_cast< int >(0xFFFFFFFF)))));
    }

    inline InaVecMaskAVX512COMMON< float > isNegativeMask() const {
        return _mm512_castps_si512(_mm512_maskz_mov_ps(_mm512_cmp_ps_mask(_mm512_setzero_ps(), vec, _CMP_GE_OQ),
                                                       _mm512_castsi512_ps(_mm512_set1_epi32(static_cast< int >(0xFFFFFFFF)))));
    }

    inline InaVecMaskAVX512COMMON< float > isPositiveStrictMask() const {
        return _mm512_castps_si512(_mm512_maskz_mov_ps(_mm512_cmp_ps_mask(_mm512_setzero_ps(), vec, _CMP_LT_OQ),
                                                       _mm512_castsi512_ps(_mm512_set1_epi32(static_cast< int >(0xFFFFFFFF)))));
    }

    inline InaVecMaskAVX512COMMON< float > isNegativeStrictMask() const {
        return _mm512_castps_si512(_mm512_maskz_mov_ps(_mm512_cmp_ps_mask(_mm512_setzero_ps(), vec, _CMP_GT_OQ),
                                                       _mm512_castsi512_ps(_mm512_set1_epi32(static_cast< int >(0xFFFFFFFF)))));
    }

    inline InaVecMaskAVX512COMMON< float > isZeroMask() const {
        return _mm512_castps_si512(_mm512_maskz_mov_ps(_mm512_cmp_ps_mask(_mm512_setzero_ps(), vec, _CMP_EQ_OQ),
                                                       _mm512_castsi512_ps(_mm512_set1_epi32(static_cast< int >(0xFFFFFFFF)))));
    }

    inline InaVecMaskAVX512COMMON< float > isNotZeroMask() const {
        return _mm512_castps_si512(_mm512_maskz_mov_ps(_mm512_cmp_ps_mask(_mm512_setzero_ps(), vec, _CMP_NEQ_OQ),
                                                       _mm512_castsi512_ps(_mm512_set1_epi32(static_cast< int >(0xFFFFFFFF)))));
    }

    // Static basic methods
    inline static InaVecAVX512COMMON GetZero() {
        return InaVecAVX512COMMON(_mm512_setzero_ps());
    }

    inline static InaVecAVX512COMMON GetOne() {
        return InaVecAVX512COMMON(_mm512_set1_ps(1));
    }

    inline static InaVecAVX512COMMON Min(const InaVecAVX512COMMON& inVec1, const InaVecAVX512COMMON& inVec2) {
        return _mm512_min_ps(inVec1.vec, inVec2.vec);
    }

    inline static InaVecAVX512COMMON Max(const InaVecAVX512COMMON& inVec1, const InaVecAVX512COMMON& inVec2) {
        return _mm512_max_ps(inVec1.vec, inVec2.vec);
    }

    inline static InaVecAVX512COMMON IsLowerOrEqual(const InaVecAVX512COMMON& inVec1, const InaVecAVX512COMMON& inVec2) {
        const __mmask16 testResult = _mm512_cmp_ps_mask(inVec1.vec, inVec2.vec, _CMP_LE_OQ);
        const __m512 ones          = _mm512_set1_ps(1);
        return _mm512_maskz_mov_ps(testResult, ones);
    }

    inline static InaVecAVX512COMMON IsLower(const InaVecAVX512COMMON& inVec1, const InaVecAVX512COMMON& inVec2) {
        const __mmask16 testResult = _mm512_cmp_ps_mask(inVec1.vec, inVec2.vec, _CMP_LT_OQ);
        const __m512 ones          = _mm512_set1_ps(1);
        return _mm512_maskz_mov_ps(testResult, ones);
    }

    inline static InaVecAVX512COMMON IsGreaterOrEqual(const InaVecAVX512COMMON& inVec1, const InaVecAVX512COMMON& inVec2) {
        const __mmask16 testResult = _mm512_cmp_ps_mask(inVec1.vec, inVec2.vec, _CMP_GE_OQ);
        const __m512 ones          = _mm512_set1_ps(1);
        return _mm512_maskz_mov_ps(testResult, ones);
    }

    inline static InaVecAVX512COMMON IsGreater(const InaVecAVX512COMMON& inVec1, const InaVecAVX512COMMON& inVec2) {
        const __mmask16 testResult = _mm512_cmp_ps_mask(inVec1.vec, inVec2.vec, _CMP_GT_OQ);
        const __m512 ones          = _mm512_set1_ps(1);
        return _mm512_maskz_mov_ps(testResult, ones);
    }

    inline static InaVecAVX512COMMON IsEqual(const InaVecAVX512COMMON& inVec1, const InaVecAVX512COMMON& inVec2) {
        const __mmask16 testResult = _mm512_cmp_ps_mask(inVec1.vec, inVec2.vec, _CMP_EQ_OQ);
        const __m512 ones          = _mm512_set1_ps(1);
        return _mm512_maskz_mov_ps(testResult, ones);
    }

    inline static InaVecAVX512COMMON IsNotEqual(const InaVecAVX512COMMON& inVec1, const InaVecAVX512COMMON& inVec2) {
        const __mmask16 testResult = _mm512_cmp_ps_mask(inVec1.vec, inVec2.vec, _CMP_NEQ_OQ);
        const __m512 ones          = _mm512_set1_ps(1);
        return _mm512_maskz_mov_ps(testResult, ones);
    }

    inline static InaVecMaskAVX512COMMON< float > IsLowerOrEqualMask(const InaVecAVX512COMMON& inVec1, const InaVecAVX512COMMON& inVec2) {
        return _mm512_castps_si512(_mm512_maskz_mov_ps(_mm512_cmp_ps_mask(inVec1.vec, inVec2.vec, _CMP_LE_OQ),
                                                       _mm512_castsi512_ps(_mm512_set1_epi32(static_cast< int >(0xFFFFFFFF)))));
    }

    inline static InaVecMaskAVX512COMMON< float > IsLowerMask(const InaVecAVX512COMMON& inVec1, const InaVecAVX512COMMON& inVec2) {
        return _mm512_castps_si512(_mm512_maskz_mov_ps(_mm512_cmp_ps_mask(inVec1.vec, inVec2.vec, _CMP_LT_OQ),
                                                       _mm512_castsi512_ps(_mm512_set1_epi32(static_cast< int >(0xFFFFFFFF)))));
    }

    inline static InaVecMaskAVX512COMMON< float > IsGreaterOrEqualMask(const InaVecAVX512COMMON& inVec1, const InaVecAVX512COMMON& inVec2) {
        return _mm512_castps_si512(_mm512_maskz_mov_ps(_mm512_cmp_ps_mask(inVec1.vec, inVec2.vec, _CMP_GE_OQ),
                                                       _mm512_castsi512_ps(_mm512_set1_epi32(static_cast< int >(0xFFFFFFFF)))));
    }

    inline static InaVecMaskAVX512COMMON< float > IsGreaterMask(const InaVecAVX512COMMON& inVec1, const InaVecAVX512COMMON& inVec2) {
        return _mm512_castps_si512(_mm512_maskz_mov_ps(_mm512_cmp_ps_mask(inVec1.vec, inVec2.vec, _CMP_GT_OQ),
                                                       _mm512_castsi512_ps(_mm512_set1_epi32(static_cast< int >(0xFFFFFFFF)))));
    }

    inline static InaVecMaskAVX512COMMON< float > IsEqualMask(const InaVecAVX512COMMON& inVec1, const InaVecAVX512COMMON& inVec2) {
        return _mm512_castps_si512(_mm512_maskz_mov_ps(_mm512_cmp_ps_mask(inVec1.vec, inVec2.vec, _CMP_EQ_OQ),
                                                       _mm512_castsi512_ps(_mm512_set1_epi32(static_cast< int >(0xFFFFFFFF)))));
    }

    inline static InaVecMaskAVX512COMMON< float > IsNotEqualMask(const InaVecAVX512COMMON& inVec1, const InaVecAVX512COMMON& inVec2) {
        return _mm512_castps_si512(_mm512_maskz_mov_ps(_mm512_cmp_ps_mask(inVec1.vec, inVec2.vec, _CMP_NEQ_OQ),
                                                       _mm512_castsi512_ps(_mm512_set1_epi32(static_cast< int >(0xFFFFFFFF)))));
    }

    inline static InaVecAVX512COMMON BitsAnd(const InaVecAVX512COMMON& inVec1, const InaVecAVX512COMMON& inVec2) {
        return _mm512_castsi512_ps(_mm512_and_si512(_mm512_castps_si512(inVec1.vec), _mm512_castps_si512(inVec2.vec)));
    }

    inline static InaVecAVX512COMMON BitsNotAnd(const InaVecAVX512COMMON& inVec1, const InaVecAVX512COMMON& inVec2) {
        return _mm512_castsi512_ps(_mm512_andnot_si512(_mm512_castps_si512(inVec1.vec), _mm512_castps_si512(inVec2.vec)));
    }

    inline static InaVecAVX512COMMON BitsOr(const InaVecAVX512COMMON& inVec1, const InaVecAVX512COMMON& inVec2) {
        return _mm512_castsi512_ps(_mm512_or_si512(_mm512_castps_si512(inVec1.vec), _mm512_castps_si512(inVec2.vec)));
    }

    inline static InaVecAVX512COMMON BitsXor(const InaVecAVX512COMMON& inVec1, const InaVecAVX512COMMON& inVec2) {
        return _mm512_castsi512_ps(_mm512_xor_si512(_mm512_castps_si512(inVec1.vec), _mm512_castps_si512(inVec2.vec)));
    }

    inline static const char* GetName() {
        return "InaVecAVX512COMMON<float>";
    }

    inline static InaIfElse< InaVecAVX512COMMON< float > >::ThenClass If(const InaVecMaskAVX512COMMON< float >& inTest) {
        return InaIfElse< InaVecAVX512COMMON< float > >::IfClass().If(inTest);
    }

    inline static InaVecAVX512COMMON IfElse(const InaVecMaskAVX512COMMON< float >& inMask, const InaVecAVX512COMMON& inIfTrue, const InaVecAVX512COMMON& inIfFalse) {
        return _mm512_castsi512_ps(_mm512_or_si512(_mm512_castps_si512(IfTrue(inMask, inIfTrue.vec).vec),
                                                   _mm512_castps_si512(IfFalse(inMask, inIfFalse.vec).vec)));
    }

    inline static InaVecAVX512COMMON IfTrue(const InaVecMaskAVX512COMMON< float >& inMask, const InaVecAVX512COMMON& inIfTrue) {
        return _mm512_castsi512_ps(_mm512_and_si512(inMask.getMask(), _mm512_castps_si512(inIfTrue.vec)));
    }

    inline static InaVecAVX512COMMON IfFalse(const InaVecMaskAVX512COMMON< float >& inMask, const InaVecAVX512COMMON& inIfFalse) {
        return _mm512_castsi512_ps(_mm512_andnot_si512(inMask.getMask(), _mm512_castps_si512(inIfFalse.vec)));
    }

    // Inner operators
    inline InaVecAVX512COMMON< float >& operator+=(const InaVecAVX512COMMON< float >& inVec) {
        vec = _mm512_add_ps(vec, inVec.vec);
        return *this;
    }

    inline InaVecAVX512COMMON< float >& operator-=(const InaVecAVX512COMMON< float >& inVec) {
        vec = _mm512_sub_ps(vec, inVec.vec);
        return *this;
    }

    inline InaVecAVX512COMMON< float >& operator/=(const InaVecAVX512COMMON< float >& inVec) {
        vec = _mm512_div_ps(vec, inVec.vec);
        return *this;
    }

    inline InaVecAVX512COMMON< float >& operator*=(const InaVecAVX512COMMON< float >& inVec) {
        vec = _mm512_mul_ps(vec, inVec.vec);
        return *this;
    }

    inline InaVecAVX512COMMON< float > operator-() const {
        const __m512 minus0 = _mm512_castsi512_ps(_mm512_set1_epi32(static_cast< int >(0x80000000)));
        return _mm512_castsi512_ps(_mm512_xor_si512(_mm512_castps_si512(vec), _mm512_castps_si512(minus0)));
    }

    inline InaVecAVX512COMMON< float > pow(std::size_t power) const {
        return InaUtils::FastPow< InaVecAVX512COMMON< float > >(*this, power);
    }
};

// Bits operators
inline InaVecAVX512COMMON< float > operator&(const InaVecAVX512COMMON< float >& inVec1, const InaVecAVX512COMMON< float >& inVec2) {
    return InaVecAVX512COMMON< float >::BitsAnd(inVec1, inVec2);
}

inline InaVecAVX512COMMON< float > operator|(const InaVecAVX512COMMON< float >& inVec1, const InaVecAVX512COMMON< float >& inVec2) {
    return InaVecAVX512COMMON< float >::BitsOr(inVec1, inVec2);
}

inline InaVecAVX512COMMON< float > operator^(const InaVecAVX512COMMON< float >& inVec1, const InaVecAVX512COMMON< float >& inVec2) {
    return InaVecAVX512COMMON< float >::BitsXor(inVec1, inVec2);
}

// Dual operators
inline InaVecAVX512COMMON< float > operator+(const InaVecAVX512COMMON< float >& inVec1, const InaVecAVX512COMMON< float >& inVec2) {
    return _mm512_add_ps(inVec1.getVec(), inVec2.getVec());
}

inline InaVecAVX512COMMON< float > operator-(const InaVecAVX512COMMON< float >& inVec1, const InaVecAVX512COMMON< float >& inVec2) {
    return _mm512_sub_ps(inVec1.getVec(), inVec2.getVec());
}

inline InaVecAVX512COMMON< float > operator/(const InaVecAVX512COMMON< float >& inVec1, const InaVecAVX512COMMON< float >& inVec2) {
    return _mm512_div_ps(inVec1.getVec(), inVec2.getVec());
}

inline InaVecAVX512COMMON< float > operator*(const InaVecAVX512COMMON< float >& inVec1, const InaVecAVX512COMMON< float >& inVec2) {
    return _mm512_mul_ps(inVec1.getVec(), inVec2.getVec());
}

// Tests and comparions
inline InaVecMaskAVX512COMMON< float > operator<(const InaVecAVX512COMMON< float >& inVec1, const InaVecAVX512COMMON< float >& inVec2) {
    return InaVecAVX512COMMON< float >::IsLowerMask(inVec1, inVec2);
}

inline InaVecMaskAVX512COMMON< float > operator<=(const InaVecAVX512COMMON< float >& inVec1, const InaVecAVX512COMMON< float >& inVec2) {
    return InaVecAVX512COMMON< float >::IsLowerOrEqualMask(inVec1, inVec2);
}

inline InaVecMaskAVX512COMMON< float > operator>(const InaVecAVX512COMMON< float >& inVec1, const InaVecAVX512COMMON< float >& inVec2) {
    return InaVecAVX512COMMON< float >::IsGreaterMask(inVec1, inVec2);
}

inline InaVecMaskAVX512COMMON< float > operator>=(const InaVecAVX512COMMON< float >& inVec1, const InaVecAVX512COMMON< float >& inVec2) {
    return InaVecAVX512COMMON< float >::IsGreaterOrEqualMask(inVec1, inVec2);
}

inline InaVecMaskAVX512COMMON< float > operator==(const InaVecAVX512COMMON< float >& inVec1, const InaVecAVX512COMMON< float >& inVec2) {
    return InaVecAVX512COMMON< float >::IsEqualMask(inVec1, inVec2);
}

inline InaVecMaskAVX512COMMON< float > operator!=(const InaVecAVX512COMMON< float >& inVec1, const InaVecAVX512COMMON< float >& inVec2) {
    return InaVecAVX512COMMON< float >::IsNotEqualMask(inVec1, inVec2);
}


#endif
