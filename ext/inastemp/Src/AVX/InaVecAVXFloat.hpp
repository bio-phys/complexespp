///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef INAVECAVXFLOAT_HPP
#define INAVECAVXFLOAT_HPP

#include "InastempConfig.h"
#include "Common/InaIfElse.hpp"
#include "Common/InaUtils.hpp"

#ifndef INASTEMP_USE_AVX
#error InaVecAVX<float> is included but AVX is not enable in the configuration
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
class alignas(32) InaVecMaskAVX< float > {
    __m256i mask;

public:
    // Classic constructors
    inline InaVecMaskAVX() {
    }

    InaVecMaskAVX(const InaVecMaskAVX&) = default;
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
        mask = (inBool ? _mm256_set1_epi32(static_cast< int >(0xFFFFFFFF)) : _mm256_setzero_si256());
    }

    inline InaVecMaskAVX& operator=(const bool inBool) {
        mask = (inBool ? _mm256_set1_epi32(static_cast< int >(0xFFFFFFFF)) : _mm256_setzero_si256());
        return (*this);
    }

    // Binary methods
    inline InaVecMaskAVX Not() const {
        return NotAnd(mask, _mm256_set1_epi32(static_cast< int >(0xFFFFFFFF)));
    }

    inline bool isAllTrue() const {
        // true if all FF => !FF => 0 & FF => 0
        return _mm256_testc_si256(mask, _mm256_set1_epi32(static_cast< int >(0xFFFFFFFF)));
    }

    inline bool isAllFalse() const {
        // true if all zero
        return _mm256_testz_si256(mask, mask);
    }

    // Double args methods
    inline static InaVecMaskAVX And(const InaVecMaskAVX& inMask1, const InaVecMaskAVX& inMask2) {
        // AVX2 return InaVecMaskAVX(_mm256_and_si256(inMask1.mask, inMask2.mask));
        return InaVecMaskAVX(_mm256_castps_si256(_mm256_and_ps(_mm256_castsi256_ps(inMask1.mask),
                                                               _mm256_castsi256_ps(inMask2.mask))));
    }

    inline static InaVecMaskAVX NotAnd(const InaVecMaskAVX& inMask1, const InaVecMaskAVX& inMask2) {
        // AVX2 return InaVecMaskAVX(_mm256_andnot_si256(inMask1.mask, inMask2.mask));
        return InaVecMaskAVX(_mm256_castps_si256(_mm256_andnot_ps(_mm256_castsi256_ps(inMask1.mask),
                                                                  _mm256_castsi256_ps(inMask2.mask))));
    }

    inline static InaVecMaskAVX Or(const InaVecMaskAVX& inMask1, const InaVecMaskAVX& inMask2) {
        // AVX2 return InaVecMaskAVX(_mm256_or_si256(inMask1.mask, inMask2.mask));
        return InaVecMaskAVX(_mm256_castps_si256(_mm256_or_ps(_mm256_castsi256_ps(inMask1.mask),
                                                              _mm256_castsi256_ps(inMask2.mask))));
    }

    inline static InaVecMaskAVX Xor(const InaVecMaskAVX& inMask1, const InaVecMaskAVX& inMask2) {
        // AVX2 return InaVecMaskAVX(_mm256_xor_si256(inMask1.mask, inMask2.mask));
        return InaVecMaskAVX(_mm256_castps_si256(_mm256_xor_ps(_mm256_castsi256_ps(inMask1.mask),
                                                               _mm256_castsi256_ps(inMask2.mask))));
    }

    inline static bool IsEqual(const InaVecMaskAVX& inMask1, const InaVecMaskAVX& inMask2) {
        return _mm256_testz_si256(_mm256_castps_si256(_mm256_xor_ps(_mm256_castsi256_ps(inMask1.mask),
                                                                    _mm256_castsi256_ps(inMask2.mask))),
                                  _mm256_set1_epi32(static_cast< int >(0xFFFFFFFF))); // return CF
    }

    inline static bool IsNotEqual(const InaVecMaskAVX& inMask1, const InaVecMaskAVX& inMask2) {
        return !_mm256_testz_si256(_mm256_castps_si256(_mm256_xor_ps(_mm256_castsi256_ps(inMask1.mask),
                                                                     _mm256_castsi256_ps(inMask2.mask))),
                                   _mm256_set1_epi32(static_cast< int >(0xFFFFFFFF))); // return CF
    }
};

// Mask must have operators
inline InaVecMaskAVX< float > operator&(const InaVecMaskAVX< float >& inMask1, const InaVecMaskAVX< float >& inMask2) {
    return InaVecMaskAVX< float >::And(inMask1, inMask2);
}

inline InaVecMaskAVX< float > operator|(const InaVecMaskAVX< float >& inMask1, const InaVecMaskAVX< float >& inMask2) {
    return InaVecMaskAVX< float >::Or(inMask1, inMask2);
}

inline InaVecMaskAVX< float > operator^(const InaVecMaskAVX< float >& inMask1, const InaVecMaskAVX< float >& inMask2) {
    return InaVecMaskAVX< float >::Xor(inMask1, inMask2);
}

inline bool operator==(const InaVecMaskAVX< float >& inMask1, const InaVecMaskAVX< float >& inMask2) {
    return InaVecMaskAVX< float >::IsEqual(inMask1, inMask2);
}

inline bool operator!=(const InaVecMaskAVX< float >& inMask1, const InaVecMaskAVX< float >& inMask2) {
    return InaVecMaskAVX< float >::IsNotEqual(inMask1, inMask2);
}

// Vec type
template <>
class alignas(32) InaVecAVX< float > {
protected:
    __m256 vec;

public:
    using VecRawType            = __m256;
    using MaskType              = InaVecMaskAVX< float >;
    using RealType              = float;
    static const int VecLength  = 8;
    static const int Alignement = 32;

    inline InaVecAVX() {
    }
    inline InaVecAVX(const InaVecAVX&) = default;
    inline InaVecAVX& operator=(const InaVecAVX&) = default;

    // Constructor from raw type
    inline /*not explicit*/ InaVecAVX(const __m256 inVec)
    : vec(inVec) {
    }

    inline InaVecAVX& operator=(const __m256 inVec) {
        vec = inVec;
        return *this;
    }

    inline void setFromRawType(const __m256 inVec) {
        vec = inVec;
    }

    inline explicit operator __m256() const {
        return vec;
    }

    inline __m256 getVec() const {
        return vec;
    }

    // Constructor from scalar
    inline /*not explicit*/ InaVecAVX(const float val)
    : vec(_mm256_set1_ps(val)) {
    }

    inline InaVecAVX& operator=(const float val) {
        vec = _mm256_set1_ps(val);
        return *this;
    }

    inline void setFromScalar(const float val) {
        vec = _mm256_set1_ps(val);
    }

    // Constructor from vec
    inline InaVecAVX(const std::initializer_list< float > lst)
    : InaVecAVX(lst.begin()) {
    }

    inline explicit InaVecAVX(const float ptr[])
    : vec(_mm256_loadu_ps(ptr)) {
    }

    inline InaVecAVX& setFromArray(const float ptr[]) {
        vec = _mm256_loadu_ps(ptr);
        return *this;
    }

    inline InaVecAVX& setFromAlignedArray(const float ptr[]) {
        vec = _mm256_load_ps(ptr);
        return *this;
    }

    inline InaVecAVX& setFromIndirectArray(const float values[], const int inIndirection[]) {
        vec = _mm256_set_ps(
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

    inline InaVecAVX& setFromIndirect2DArray(const float inArray[], const int inIndirection1[],
                                             const int inLeadingDimension, const int inIndirection2[]) {
        vec = _mm256_set_ps(
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
        _mm256_storeu_ps(ptr, vec);
    }

    inline void storeInAlignedArray(float ptr[]) const {
        _mm256_store_ps(ptr, vec);
    }

    // Acce to individual values
    inline float at(const int index) const {
        alignas(Alignement) float allval[VecLength];
        _mm256_store_ps(allval, vec);
        return allval[index];
    }

    // Horizontal operation
    inline float horizontalSum() const {
        const __m128 valupper = _mm256_extractf128_ps(vec, 1);
        const __m128 rest     = _mm256_extractf128_ps(vec, 0);
        const __m128 valval   = _mm_add_ps(valupper,
                                         rest);
        __m128 valsum = _mm_add_ps(_mm_permute_ps(valval, 0x1B), valval);
        __m128 res    = _mm_add_ps(_mm_permute_ps(valsum, 0xB1), valsum);
        return _mm_cvtss_f32(res);
    }

    inline float horizontalMul() const {
        const __m128 valupper = _mm256_extractf128_ps(vec, 1);
        const __m128 rest     = _mm256_extractf128_ps(vec, 0);
        const __m128 valval   = _mm_mul_ps(valupper,
                                         rest);
        __m128 valsum = _mm_mul_ps(_mm_permute_ps(valval, 0x1B), valval);
        __m128 res    = _mm_mul_ps(_mm_permute_ps(valsum, 0xB1), valsum);
        return _mm_cvtss_f32(res);
    }

    inline InaVecAVX sqrt() const {
        return _mm256_sqrt_ps(vec);
    }

    inline InaVecAVX exp() const {
#ifdef __INTEL_COMPILER
        return _mm256_exp_ps(vec);
#else
        const __m256 COEFF_LOG2E = _mm256_set1_ps(float(InaFastExp::CoeffLog2E()));
        const __m256 COEFF_A     = _mm256_set1_ps(float(InaFastExp::CoeffA32()));
        const __m256 COEFF_B     = _mm256_set1_ps(float(InaFastExp::CoeffB32()));
        const __m256 COEFF_P5_A  = _mm256_set1_ps(float(InaFastExp::GetCoefficient6_5()));
        const __m256 COEFF_P5_B  = _mm256_set1_ps(float(InaFastExp::GetCoefficient6_4()));
        const __m256 COEFF_P5_C  = _mm256_set1_ps(float(InaFastExp::GetCoefficient6_3()));
        const __m256 COEFF_P5_D  = _mm256_set1_ps(float(InaFastExp::GetCoefficient6_2()));
        const __m256 COEFF_P5_E  = _mm256_set1_ps(float(InaFastExp::GetCoefficient6_1()));
        const __m256 COEFF_P5_F  = _mm256_set1_ps(float(InaFastExp::GetCoefficient6_0()));

        __m256 x = _mm256_mul_ps(vec, COEFF_LOG2E);

        const __m256 fractional_part = _mm256_sub_ps(x, InaVecAVX(x).floor().vec);

        __m256 factor = _mm256_add_ps(_mm256_mul_ps(_mm256_add_ps(_mm256_mul_ps(_mm256_add_ps(
                                                                                    _mm256_mul_ps(_mm256_add_ps(_mm256_mul_ps(_mm256_add_ps(_mm256_mul_ps(
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

        x = _mm256_sub_ps(x, factor);

        __m256i castedInteger = _mm256_cvtps_epi32(_mm256_add_ps(_mm256_mul_ps(COEFF_A, x), COEFF_B));

        return _mm256_castsi256_ps(castedInteger);
#endif
    }

    inline InaVecAVX expLowAcc() const {
        const __m256 COEFF_LOG2E = _mm256_set1_ps(float(InaFastExp::CoeffLog2E()));
        const __m256 COEFF_A     = _mm256_set1_ps(float(InaFastExp::CoeffA32()));
        const __m256 COEFF_B     = _mm256_set1_ps(float(InaFastExp::CoeffB32()));
        const __m256 COEFF_P5_D  = _mm256_set1_ps(float(InaFastExp::GetCoefficient3_2()));
        const __m256 COEFF_P5_E  = _mm256_set1_ps(float(InaFastExp::GetCoefficient3_1()));
        const __m256 COEFF_P5_F  = _mm256_set1_ps(float(InaFastExp::GetCoefficient3_0()));

        __m256 x = _mm256_mul_ps(vec, COEFF_LOG2E);

        const __m256 fractional_part = _mm256_sub_ps(x, InaVecAVX(x).floor().vec);

        __m256 factor = _mm256_add_ps(_mm256_mul_ps(
                                          _mm256_add_ps(_mm256_mul_ps(
                                                            COEFF_P5_D, fractional_part),
                                                        COEFF_P5_E),
                                          fractional_part),
                                      COEFF_P5_F);

        x = _mm256_sub_ps(x, factor);

        __m256i castedInteger = _mm256_cvtps_epi32(_mm256_add_ps(_mm256_mul_ps(COEFF_A, x), COEFF_B));

        return _mm256_castsi256_ps(castedInteger);
    }

    inline InaVecAVX rsqrt() const {
        return _mm256_set1_ps(1) / _mm256_sqrt_ps(vec); // _mm256_rsqrt_ps(val); not accurate enough
    }

    inline InaVecAVX abs() const {
        const __m256 minus0 = _mm256_castsi256_ps(_mm256_set1_epi32(static_cast< int >(0x80000000)));
        return _mm256_andnot_ps(minus0, vec);
    }

    inline InaVecAVX floor() const {
        return _mm256_floor_ps(vec);
    }

    inline InaVecAVX signOf() const {
        const __m256 minus0 = _mm256_castsi256_ps(_mm256_set1_epi32(static_cast< int >(0x80000000)));
        const __m256 signs  = _mm256_and_ps(vec, minus0);
        return _mm256_and_ps(_mm256_cmp_ps(_mm256_setzero_ps(), vec, _CMP_NEQ_OQ),
                             _mm256_or_ps(signs, _mm256_set1_ps(1)));
    }

    inline InaVecAVX isPositive() const {
        const __m256 greater = _mm256_cmp_ps(_mm256_setzero_ps(), vec, _CMP_LE_OQ);
        const __m256 ones    = _mm256_set1_ps(1);
        return _mm256_and_ps(greater, ones);
    }

    inline InaVecAVX isNegative() const {
        const __m256 less = _mm256_cmp_ps(_mm256_setzero_ps(), vec, _CMP_GE_OQ);
        const __m256 ones = _mm256_set1_ps(1);
        return _mm256_and_ps(less, ones);
    }

    inline InaVecAVX isPositiveStrict() const {
        const __m256 greater = _mm256_cmp_ps(_mm256_setzero_ps(), vec, _CMP_LT_OQ);
        const __m256 ones    = _mm256_set1_ps(1);
        return _mm256_and_ps(greater, ones);
    }

    inline InaVecAVX isNegativeStrict() const {
        const __m256 less = _mm256_cmp_ps(_mm256_setzero_ps(), vec, _CMP_GT_OQ);
        const __m256 ones = _mm256_set1_ps(1);
        return _mm256_and_ps(less, ones);
    }

    inline InaVecAVX isZero() const {
        const __m256 equalZero = _mm256_cmp_ps(_mm256_setzero_ps(), vec, _CMP_EQ_OQ);
        const __m256 ones      = _mm256_set1_ps(1);
        return _mm256_and_ps(equalZero, ones);
    }

    inline InaVecAVX isNotZero() const {
        const __m256 equalZero = _mm256_cmp_ps(_mm256_setzero_ps(), vec, _CMP_NEQ_OQ);
        const __m256 ones      = _mm256_set1_ps(1);
        return _mm256_and_ps(equalZero, ones);
    }

    inline InaVecMaskAVX< float > isPositiveMask() const {
        return _mm256_castps_si256(_mm256_cmp_ps(_mm256_setzero_ps(), vec, _CMP_LE_OQ));
    }

    inline InaVecMaskAVX< float > isNegativeMask() const {
        return _mm256_castps_si256(_mm256_cmp_ps(_mm256_setzero_ps(), vec, _CMP_GE_OQ));
    }

    inline InaVecMaskAVX< float > isPositiveStrictMask() const {
        return _mm256_castps_si256(_mm256_cmp_ps(_mm256_setzero_ps(), vec, _CMP_LT_OQ));
    }

    inline InaVecMaskAVX< float > isNegativeStrictMask() const {
        return _mm256_castps_si256(_mm256_cmp_ps(_mm256_setzero_ps(), vec, _CMP_GT_OQ));
    }

    inline InaVecMaskAVX< float > isZeroMask() const {
        return _mm256_castps_si256(_mm256_cmp_ps(_mm256_setzero_ps(), vec, _CMP_EQ_OQ));
    }

    inline InaVecMaskAVX< float > isNotZeroMask() const {
        return _mm256_castps_si256(_mm256_cmp_ps(_mm256_setzero_ps(), vec, _CMP_NEQ_OQ));
    }

    // Static basic methods
    inline static InaVecAVX GetZero() {
        return InaVecAVX(_mm256_setzero_ps());
    }

    inline static InaVecAVX GetOne() {
        return InaVecAVX(_mm256_set1_ps(1));
    }

    inline static InaVecAVX Min(const InaVecAVX& inVec1, const InaVecAVX& inVec2) {
        return _mm256_min_ps(inVec1.vec, inVec2.vec);
    }

    inline static InaVecAVX Max(const InaVecAVX& inVec1, const InaVecAVX& inVec2) {
        return _mm256_max_ps(inVec1.vec, inVec2.vec);
    }

    inline static InaVecAVX IsLowerOrEqual(const InaVecAVX& inVec1, const InaVecAVX& inVec2) {
        const __m256 testResult = _mm256_cmp_ps(inVec1.vec, inVec2.vec, _CMP_LE_OQ);
        const __m256 ones       = _mm256_set1_ps(1);
        return _mm256_and_ps(testResult, ones);
    }

    inline static InaVecAVX IsLower(const InaVecAVX& inVec1, const InaVecAVX& inVec2) {
        const __m256 testResult = _mm256_cmp_ps(inVec1.vec, inVec2.vec, _CMP_LT_OQ);
        const __m256 ones       = _mm256_set1_ps(1);
        return _mm256_and_ps(testResult, ones);
    }

    inline static InaVecAVX IsGreaterOrEqual(const InaVecAVX& inVec1, const InaVecAVX& inVec2) {
        const __m256 testResult = _mm256_cmp_ps(inVec1.vec, inVec2.vec, _CMP_GE_OQ);
        const __m256 ones       = _mm256_set1_ps(1);
        return _mm256_and_ps(testResult, ones);
    }

    inline static InaVecAVX IsGreater(const InaVecAVX& inVec1, const InaVecAVX& inVec2) {
        const __m256 testResult = _mm256_cmp_ps(inVec1.vec, inVec2.vec, _CMP_GT_OQ);
        const __m256 ones       = _mm256_set1_ps(1);
        return _mm256_and_ps(testResult, ones);
    }

    inline static InaVecAVX IsEqual(const InaVecAVX& inVec1, const InaVecAVX& inVec2) {
        const __m256 testResult = _mm256_cmp_ps(inVec1.vec, inVec2.vec, _CMP_EQ_OQ);
        const __m256 ones       = _mm256_set1_ps(1);
        return _mm256_and_ps(testResult, ones);
    }

    inline static InaVecAVX IsNotEqual(const InaVecAVX& inVec1, const InaVecAVX& inVec2) {
        const __m256 testResult = _mm256_cmp_ps(inVec1.vec, inVec2.vec, _CMP_NEQ_OQ);
        const __m256 ones       = _mm256_set1_ps(1);
        return _mm256_and_ps(testResult, ones);
    }

    inline static InaVecMaskAVX< float > IsLowerOrEqualMask(const InaVecAVX& inVec1, const InaVecAVX& inVec2) {
        return _mm256_castps_si256(_mm256_cmp_ps(inVec1.vec, inVec2.vec, _CMP_LE_OQ));
    }

    inline static InaVecMaskAVX< float > IsLowerMask(const InaVecAVX& inVec1, const InaVecAVX& inVec2) {
        return _mm256_castps_si256(_mm256_cmp_ps(inVec1.vec, inVec2.vec, _CMP_LT_OQ));
    }

    inline static InaVecMaskAVX< float > IsGreaterOrEqualMask(const InaVecAVX& inVec1, const InaVecAVX& inVec2) {
        return _mm256_castps_si256(_mm256_cmp_ps(inVec1.vec, inVec2.vec, _CMP_GE_OQ));
    }

    inline static InaVecMaskAVX< float > IsGreaterMask(const InaVecAVX& inVec1, const InaVecAVX& inVec2) {
        return _mm256_castps_si256(_mm256_cmp_ps(inVec1.vec, inVec2.vec, _CMP_GT_OQ));
    }

    inline static InaVecMaskAVX< float > IsEqualMask(const InaVecAVX& inVec1, const InaVecAVX& inVec2) {
        return _mm256_castps_si256(_mm256_cmp_ps(inVec1.vec, inVec2.vec, _CMP_EQ_OQ));
    }

    inline static InaVecMaskAVX< float > IsNotEqualMask(const InaVecAVX& inVec1, const InaVecAVX& inVec2) {
        return _mm256_castps_si256(_mm256_cmp_ps(inVec1.vec, inVec2.vec, _CMP_NEQ_OQ));
    }

    inline static InaVecAVX BitsAnd(const InaVecAVX& inVec1, const InaVecAVX& inVec2) {
        return _mm256_and_ps(inVec1.vec, inVec2.vec);
    }

    inline static InaVecAVX BitsNotAnd(const InaVecAVX& inVec1, const InaVecAVX& inVec2) {
        return _mm256_andnot_ps(inVec1.vec, inVec2.vec);
    }

    inline static InaVecAVX BitsOr(const InaVecAVX& inVec1, const InaVecAVX& inVec2) {
        return _mm256_or_ps(inVec1.vec, inVec2.vec);
    }

    inline static InaVecAVX BitsXor(const InaVecAVX& inVec1, const InaVecAVX& inVec2) {
        return _mm256_xor_ps(inVec1.vec, inVec2.vec);
    }

    inline static const char* GetName() {
        return "InaVecAVX<float>";
    }

    inline static InaIfElse< InaVecAVX< float > >::ThenClass If(const InaVecMaskAVX< float >& inTest) {
        return InaIfElse< InaVecAVX< float > >::IfClass().If(inTest);
    }

    inline static InaVecAVX IfElse(const InaVecMaskAVX< float >& inMask, const InaVecAVX& inIfTrue, const InaVecAVX& inIfFalse) {
        return _mm256_or_ps(IfTrue(inMask, inIfTrue.vec).vec,
                            IfFalse(inMask, inIfFalse.vec).vec);
    }

    inline static InaVecAVX IfTrue(const InaVecMaskAVX< float >& inMask, const InaVecAVX& inIfTrue) {
        return _mm256_and_ps(_mm256_castsi256_ps(inMask.getMask()), inIfTrue.vec);
    }

    inline static InaVecAVX IfFalse(const InaVecMaskAVX< float >& inMask, const InaVecAVX& inIfFalse) {
        return _mm256_andnot_ps(_mm256_castsi256_ps(inMask.getMask()), inIfFalse.vec);
    }

    // Inner operators
    inline InaVecAVX< float >& operator+=(const InaVecAVX< float >& inVec) {
        vec = _mm256_add_ps(vec, inVec.vec);
        return *this;
    }

    inline InaVecAVX< float >& operator-=(const InaVecAVX< float >& inVec) {
        vec = _mm256_sub_ps(vec, inVec.vec);
        return *this;
    }

    inline InaVecAVX< float >& operator/=(const InaVecAVX< float >& inVec) {
        vec = _mm256_div_ps(vec, inVec.vec);
        return *this;
    }

    inline InaVecAVX< float >& operator*=(const InaVecAVX< float >& inVec) {
        vec = _mm256_mul_ps(vec, inVec.vec);
        return *this;
    }

    inline InaVecAVX< float > operator-() const {
        const __m256 minus0 = _mm256_castsi256_ps(_mm256_set1_epi32(static_cast< int >(0x80000000)));
        return _mm256_xor_ps(vec, minus0);
    }

    inline InaVecAVX< float > pow(std::size_t power) const {
        return InaUtils::FastPow< InaVecAVX< float > >(*this, power);
    }
};

// Bits operators
inline InaVecAVX< float > operator&(const InaVecAVX< float >& inVec1, const InaVecAVX< float >& inVec2) {
    return InaVecAVX< float >::BitsAnd(inVec1, inVec2);
}

inline InaVecAVX< float > operator|(const InaVecAVX< float >& inVec1, const InaVecAVX< float >& inVec2) {
    return InaVecAVX< float >::BitsOr(inVec1, inVec2);
}

inline InaVecAVX< float > operator^(const InaVecAVX< float >& inVec1, const InaVecAVX< float >& inVec2) {
    return InaVecAVX< float >::BitsXor(inVec1, inVec2);
}

// Dual operators
inline InaVecAVX< float > operator+(const InaVecAVX< float >& inVec1, const InaVecAVX< float >& inVec2) {
    return _mm256_add_ps(inVec1.getVec(), inVec2.getVec());
}

inline InaVecAVX< float > operator-(const InaVecAVX< float >& inVec1, const InaVecAVX< float >& inVec2) {
    return _mm256_sub_ps(inVec1.getVec(), inVec2.getVec());
}

inline InaVecAVX< float > operator/(const InaVecAVX< float >& inVec1, const InaVecAVX< float >& inVec2) {
    return _mm256_div_ps(inVec1.getVec(), inVec2.getVec());
}

inline InaVecAVX< float > operator*(const InaVecAVX< float >& inVec1, const InaVecAVX< float >& inVec2) {
    return _mm256_mul_ps(inVec1.getVec(), inVec2.getVec());
}

// Tests and comparions
inline InaVecMaskAVX< float > operator<(const InaVecAVX< float >& inVec1, const InaVecAVX< float >& inVec2) {
    return InaVecAVX< float >::IsLowerMask(inVec1, inVec2);
}

inline InaVecMaskAVX< float > operator<=(const InaVecAVX< float >& inVec1, const InaVecAVX< float >& inVec2) {
    return InaVecAVX< float >::IsLowerOrEqualMask(inVec1, inVec2);
}

inline InaVecMaskAVX< float > operator>(const InaVecAVX< float >& inVec1, const InaVecAVX< float >& inVec2) {
    return InaVecAVX< float >::IsGreaterMask(inVec1, inVec2);
}

inline InaVecMaskAVX< float > operator>=(const InaVecAVX< float >& inVec1, const InaVecAVX< float >& inVec2) {
    return InaVecAVX< float >::IsGreaterOrEqualMask(inVec1, inVec2);
}

inline InaVecMaskAVX< float > operator==(const InaVecAVX< float >& inVec1, const InaVecAVX< float >& inVec2) {
    return InaVecAVX< float >::IsEqualMask(inVec1, inVec2);
}

inline InaVecMaskAVX< float > operator!=(const InaVecAVX< float >& inVec1, const InaVecAVX< float >& inVec2) {
    return InaVecAVX< float >::IsNotEqualMask(inVec1, inVec2);
}


#endif
