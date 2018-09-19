///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef INAVECSSE3FLOAT_HPP
#define INAVECSSE3FLOAT_HPP

#include "InastempConfig.h"
#include "Common/InaIfElse.hpp"
#include "Common/InaUtils.hpp"

#ifndef INASTEMP_USE_SSE3
#error InaVecSSE3<float> is included but SSE3 is not enable in the configuration
#endif

#include "Common/InaFastExp.hpp"

#include <emmintrin.h>
#include <cmath>
#include <initializer_list>

// Forward declarations
template < class RealType >
class InaVecMaskSSE3;

template < class RealType >
class InaVecSSE3;


// Mask type
template <>
class InaVecMaskSSE3< float > {
    __m128i mask;

public:
    // Classic constructors
    inline InaVecMaskSSE3() {
    }

    inline InaVecMaskSSE3(const InaVecMaskSSE3&) = default;
    inline InaVecMaskSSE3& operator=(const InaVecMaskSSE3&) = default;

    // Native data type compatibility
    inline /*not explicit*/ InaVecMaskSSE3(const __m128i inMask)
    : mask(inMask) {
    }

    inline InaVecMaskSSE3& operator=(const __m128i inMask) {
        mask = inMask;
        return (*this);
    }

    inline explicit operator __m128i() const {
        return mask;
    }

    inline __m128i getMask() const {
        return mask;
    }

    // Bool data type compatibility
    inline explicit InaVecMaskSSE3(const bool inBool) {
        mask = (inBool ? _mm_set1_epi32(static_cast< int >(0xFFFFFFFF)) : _mm_setzero_si128());
    }

    inline InaVecMaskSSE3& operator=(const bool inBool) {
        mask = (inBool ? _mm_set1_epi32(static_cast< int >(0xFFFFFFFF)) : _mm_setzero_si128());
        return (*this);
    }

    // Binary methods
    inline InaVecMaskSSE3 Not() const {
        return NotAnd(mask, _mm_set1_epi32(static_cast< int >(0xFFFFFFFF)));
    }

    inline bool isAllTrue() const {
        // true if all FF
        return _mm_movemask_epi8(_mm_cmpeq_epi32(mask, _mm_set1_epi32(static_cast< int >(0xFFFFFFFF)))) == 0xFFFF;
    }

    inline bool isAllFalse() const {
        // true if all zero
        return _mm_movemask_epi8(_mm_cmpeq_epi32(mask, _mm_setzero_si128())) == 0xFFFF;
    }

    // Double args methods
    inline static InaVecMaskSSE3 And(const InaVecMaskSSE3& inMask1, const InaVecMaskSSE3& inMask2) {
        return InaVecMaskSSE3(_mm_and_si128(inMask1.mask, inMask2.mask));
    }

    inline static InaVecMaskSSE3 NotAnd(const InaVecMaskSSE3& inMask1, const InaVecMaskSSE3& inMask2) {
        return InaVecMaskSSE3(_mm_andnot_si128(inMask1.mask, inMask2.mask));
    }

    inline static InaVecMaskSSE3 Or(const InaVecMaskSSE3& inMask1, const InaVecMaskSSE3& inMask2) {
        return InaVecMaskSSE3(_mm_or_si128(inMask1.mask, inMask2.mask));
    }

    inline static InaVecMaskSSE3 Xor(const InaVecMaskSSE3& inMask1, const InaVecMaskSSE3& inMask2) {
        return InaVecMaskSSE3(_mm_xor_si128(inMask1.mask, inMask2.mask));
    }

    inline static bool IsEqual(const InaVecMaskSSE3& inMask1, const InaVecMaskSSE3& inMask2) {
        return _mm_movemask_epi8(_mm_cmpeq_epi32(inMask1.mask, inMask2.mask)) == 0xFFFF;
    }

    inline static bool IsNotEqual(const InaVecMaskSSE3& inMask1, const InaVecMaskSSE3& inMask2) {
        return _mm_movemask_epi8(_mm_cmpeq_epi32(inMask1.mask, inMask2.mask)) != 0xFFFF;
    }
};

// Mask must have operators
inline InaVecMaskSSE3< float > operator&(const InaVecMaskSSE3< float >& inMask1, const InaVecMaskSSE3< float >& inMask2) {
    return InaVecMaskSSE3< float >::And(inMask1, inMask2);
}

inline InaVecMaskSSE3< float > operator|(const InaVecMaskSSE3< float >& inMask1, const InaVecMaskSSE3< float >& inMask2) {
    return InaVecMaskSSE3< float >::Or(inMask1, inMask2);
}

inline InaVecMaskSSE3< float > operator^(const InaVecMaskSSE3< float >& inMask1, const InaVecMaskSSE3< float >& inMask2) {
    return InaVecMaskSSE3< float >::Xor(inMask1, inMask2);
}

inline bool operator==(const InaVecMaskSSE3< float >& inMask1, const InaVecMaskSSE3< float >& inMask2) {
    return InaVecMaskSSE3< float >::IsEqual(inMask1, inMask2);
}

inline bool operator!=(const InaVecMaskSSE3< float >& inMask1, const InaVecMaskSSE3< float >& inMask2) {
    return InaVecMaskSSE3< float >::IsNotEqual(inMask1, inMask2);
}

// Vec type
template <>
class InaVecSSE3< float > {
protected:
    __m128 vec;

public:
    using VecRawType            = __m128;
    using MaskType              = InaVecMaskSSE3< float >;
    using RealType              = float;
    static const int VecLength  = 4;
    static const int Alignement = 16;

    inline InaVecSSE3() {
    }
    inline InaVecSSE3(const InaVecSSE3&) = default;
    inline InaVecSSE3& operator=(const InaVecSSE3&) = default;

    // Constructor from raw type
    /*not explicit*/ InaVecSSE3(const __m128 inVec)
    : vec(inVec) {
    }

    inline InaVecSSE3& operator=(const __m128 inVec) {
        vec = inVec;
        return *this;
    }

    inline void setFromRawType(const __m128 inVec) {
        vec = inVec;
    }

    inline explicit operator __m128() const {
        return vec;
    }

    inline __m128 getVec() const {
        return vec;
    }

    // Constructor from scalar
    /*not explicit*/ InaVecSSE3(const float val)
    : vec(_mm_set1_ps(val)) {
    }

    inline InaVecSSE3& operator=(const float val) {
        vec = _mm_set1_ps(val);
        return *this;
    }

    inline void setFromScalar(const float val) {
        vec = _mm_set1_ps(val);
    }

    // Constructor from vec
    inline InaVecSSE3(const std::initializer_list< float > lst)
    : InaVecSSE3(lst.begin()) {
    }

    inline explicit InaVecSSE3(const float ptr[])
    : vec(_mm_loadu_ps(ptr)) {
    }

    inline InaVecSSE3& setFromArray(const float ptr[]) {
        vec = _mm_loadu_ps(ptr);
        return *this;
    }

    inline InaVecSSE3& setFromAlignedArray(const float ptr[]) {
        vec = _mm_load_ps(ptr);
        return *this;
    }

    inline InaVecSSE3& setFromIndirectArray(const float values[], const int inIndirection[]) {
        vec = _mm_set_ps(
            values[inIndirection[3]],
            values[inIndirection[2]],
            values[inIndirection[1]],
            values[inIndirection[0]]);
        return *this;
    }

    inline InaVecSSE3& setFromIndirect2DArray(const float inArray[], const int inIndirection1[],
                                              const int inLeadingDimension, const int inIndirection2[]) {
        vec = _mm_set_ps(
            inArray[inIndirection1[3] * inLeadingDimension + inIndirection2[3]],
            inArray[inIndirection1[2] * inLeadingDimension + inIndirection2[2]],
            inArray[inIndirection1[1] * inLeadingDimension + inIndirection2[1]],
            inArray[inIndirection1[0] * inLeadingDimension + inIndirection2[0]]);
        return *this;
    }

    // Move back to array
    inline void storeInArray(float ptr[]) const {
        _mm_storeu_ps(ptr, vec);
    }

    inline void storeInAlignedArray(float ptr[]) const {
        _mm_store_ps(ptr, vec);
    }

    // Acce to individual values
    inline float at(const int index) const {
        alignas(Alignement) float allval[VecLength];
        _mm_store_ps(allval, vec);
        return allval[index];
    }

    // Horizontal operation
    inline float horizontalSum() const {
        // Better than _mm_hadd_ps
        const __m128 val02_13_20_31 = _mm_add_ps(vec, _mm_movehl_ps(vec, vec));
        const __m128 res            = _mm_add_ss(val02_13_20_31, _mm_shuffle_ps(val02_13_20_31, val02_13_20_31, 1));
        return _mm_cvtss_f32(res);
    }

    inline float horizontalMul() const {
        const __m128 val02_13_20_31 = _mm_mul_ps(vec, _mm_movehl_ps(vec, vec));
        const __m128 res            = _mm_mul_ss(val02_13_20_31, _mm_shuffle_ps(val02_13_20_31, val02_13_20_31, 1));
        return _mm_cvtss_f32(res);
    }

    inline InaVecSSE3 sqrt() const {
        return _mm_sqrt_ps(vec);
    }

    inline InaVecSSE3 exp() const {
#ifdef __INTEL_COMPILER
        return _mm_exp_ps(vec);
#else
        const __m128 COEFF_LOG2E = _mm_set1_ps(float(InaFastExp::CoeffLog2E()));
        const __m128 COEFF_A     = _mm_set1_ps(float(InaFastExp::CoeffA32()));
        const __m128 COEFF_B     = _mm_set1_ps(float(InaFastExp::CoeffB32()));
        const __m128 COEFF_P5_A  = _mm_set1_ps(float(InaFastExp::GetCoefficient6_5()));
        const __m128 COEFF_P5_B  = _mm_set1_ps(float(InaFastExp::GetCoefficient6_4()));
        const __m128 COEFF_P5_C  = _mm_set1_ps(float(InaFastExp::GetCoefficient6_3()));
        const __m128 COEFF_P5_D  = _mm_set1_ps(float(InaFastExp::GetCoefficient6_2()));
        const __m128 COEFF_P5_E  = _mm_set1_ps(float(InaFastExp::GetCoefficient6_1()));
        const __m128 COEFF_P5_F  = _mm_set1_ps(float(InaFastExp::GetCoefficient6_0()));

        __m128 x = _mm_mul_ps(vec, COEFF_LOG2E);

        const __m128 fractional_part = _mm_sub_ps(x, InaVecSSE3(x).floor().vec);

        __m128 factor = _mm_add_ps(_mm_mul_ps(_mm_add_ps(_mm_mul_ps(_mm_add_ps(
                                                                        _mm_mul_ps(_mm_add_ps(_mm_mul_ps(_mm_add_ps(_mm_mul_ps(
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

        x = _mm_sub_ps(x, factor);

        __m128i castedInteger = _mm_cvtps_epi32(_mm_add_ps(_mm_mul_ps(COEFF_A, x), COEFF_B));

        return _mm_castsi128_ps(castedInteger);
#endif
    }

    inline InaVecSSE3 expLowAcc() const {
        const __m128 COEFF_LOG2E = _mm_set1_ps(float(InaFastExp::CoeffLog2E()));
        const __m128 COEFF_A     = _mm_set1_ps(float(InaFastExp::CoeffA32()));
        const __m128 COEFF_B     = _mm_set1_ps(float(InaFastExp::CoeffB32()));
        const __m128 COEFF_P5_D  = _mm_set1_ps(float(InaFastExp::GetCoefficient3_2()));
        const __m128 COEFF_P5_E  = _mm_set1_ps(float(InaFastExp::GetCoefficient3_1()));
        const __m128 COEFF_P5_F  = _mm_set1_ps(float(InaFastExp::GetCoefficient3_0()));

        __m128 x = _mm_mul_ps(vec, COEFF_LOG2E);

        const __m128 fractional_part = _mm_sub_ps(x, InaVecSSE3(x).floor().vec);

        __m128 factor = _mm_add_ps(_mm_mul_ps(
                                       _mm_add_ps(_mm_mul_ps(
                                                      COEFF_P5_D, fractional_part),
                                                  COEFF_P5_E),
                                       fractional_part),
                                   COEFF_P5_F);

        x = _mm_sub_ps(x, factor);

        __m128i castedInteger = _mm_cvtps_epi32(_mm_add_ps(_mm_mul_ps(COEFF_A, x), COEFF_B));

        return _mm_castsi128_ps(castedInteger);
    }

    inline InaVecSSE3 rsqrt() const {
        return _mm_set1_ps(1) / _mm_sqrt_ps(vec); // _mm_rsqrt_ps(val); not enough accurate
    }

    inline InaVecSSE3 abs() const {
        const __m128 minus0 = _mm_castsi128_ps(_mm_set1_epi32(static_cast< int >(0x80000000)));
        return _mm_andnot_ps(minus0, vec);
    }

    inline InaVecSSE3 floor() const {
        alignas(Alignement) float allval[VecLength];
        _mm_store_ps(allval, vec);
        for(int idx = 0; idx < VecLength; ++idx) {
            allval[idx] = std::floor(allval[idx]);
        }
        return _mm_loadu_ps(allval);
    }

    inline InaVecSSE3 signOf() const {
        const __m128 minus0 = _mm_castsi128_ps(_mm_set1_epi32(static_cast< int >(0x80000000)));
        const __m128 signs  = _mm_and_ps(vec, minus0);
        return _mm_andnot_ps(_mm_cmpeq_ps(_mm_setzero_ps(), vec), _mm_or_ps(signs, _mm_set1_ps(1)));
    }

    inline InaVecSSE3 isPositive() const {
        const __m128 greater = _mm_cmple_ps(_mm_setzero_ps(), vec);
        const __m128 ones    = _mm_set1_ps(1);
        return _mm_and_ps(greater, ones);
    }

    inline InaVecSSE3 isNegative() const {
        const __m128 less = _mm_cmpge_ps(_mm_setzero_ps(), vec);
        const __m128 ones = _mm_set1_ps(1);
        return _mm_and_ps(less, ones);
    }

    inline InaVecSSE3 isPositiveStrict() const {
        const __m128 greater = _mm_cmplt_ps(_mm_setzero_ps(), vec);
        const __m128 ones    = _mm_set1_ps(1);
        return _mm_and_ps(greater, ones);
    }

    inline InaVecSSE3 isNegativeStrict() const {
        const __m128 less = _mm_cmpgt_ps(_mm_setzero_ps(), vec);
        const __m128 ones = _mm_set1_ps(1);
        return _mm_and_ps(less, ones);
    }

    inline InaVecSSE3 isZero() const {
        const __m128 equalZero = _mm_cmpeq_ps(_mm_setzero_ps(), vec);
        const __m128 ones      = _mm_set1_ps(1);
        return _mm_and_ps(equalZero, ones);
    }

    inline InaVecSSE3 isNotZero() const {
        const __m128 equalZero = _mm_cmpeq_ps(_mm_setzero_ps(), vec);
        const __m128 ones      = _mm_set1_ps(1);
        return _mm_andnot_ps(equalZero, ones);
    }

    inline InaVecMaskSSE3< float > isPositiveMask() const {
        return _mm_castps_si128(_mm_cmple_ps(_mm_setzero_ps(), vec));
    }

    inline InaVecMaskSSE3< float > isNegativeMask() const {
        return _mm_castps_si128(_mm_cmpge_ps(_mm_setzero_ps(), vec));
    }

    inline InaVecMaskSSE3< float > isPositiveStrictMask() const {
        return _mm_castps_si128(_mm_cmplt_ps(_mm_setzero_ps(), vec));
    }

    inline InaVecMaskSSE3< float > isNegativeStrictMask() const {
        return _mm_castps_si128(_mm_cmpgt_ps(_mm_setzero_ps(), vec));
    }

    inline InaVecMaskSSE3< float > isZeroMask() const {
        return _mm_castps_si128(_mm_cmpeq_ps(_mm_setzero_ps(), vec));
    }

    inline InaVecMaskSSE3< float > isNotZeroMask() const {
        return _mm_castps_si128(_mm_cmpneq_ps(_mm_setzero_ps(), vec));
    }

    // Static basic methods
    inline static InaVecSSE3 GetZero() {
        return InaVecSSE3(_mm_setzero_ps());
    }

    inline static InaVecSSE3 GetOne() {
        return InaVecSSE3(_mm_set1_ps(1));
    }

    inline static InaVecSSE3 Min(const InaVecSSE3& inVec1, const InaVecSSE3& inVec2) {
        return _mm_min_ps(inVec1.vec, inVec2.vec);
    }

    inline static InaVecSSE3 Max(const InaVecSSE3& inVec1, const InaVecSSE3& inVec2) {
        return _mm_max_ps(inVec1.vec, inVec2.vec);
    }

    inline static InaVecSSE3 IsLowerOrEqual(const InaVecSSE3& inVec1, const InaVecSSE3& inVec2) {
        const __m128 testResult = _mm_cmple_ps(inVec1.vec, inVec2.vec);
        const __m128 ones       = _mm_set1_ps(1);
        return _mm_and_ps(testResult, ones);
    }

    inline static InaVecSSE3 IsLower(const InaVecSSE3& inVec1, const InaVecSSE3& inVec2) {
        const __m128 testResult = _mm_cmplt_ps(inVec1.vec, inVec2.vec);
        const __m128 ones       = _mm_set1_ps(1);
        return _mm_and_ps(testResult, ones);
    }

    inline static InaVecSSE3 IsGreaterOrEqual(const InaVecSSE3& inVec1, const InaVecSSE3& inVec2) {
        const __m128 testResult = _mm_cmpge_ps(inVec1.vec, inVec2.vec);
        const __m128 ones       = _mm_set1_ps(1);
        return _mm_and_ps(testResult, ones);
    }

    inline static InaVecSSE3 IsGreater(const InaVecSSE3& inVec1, const InaVecSSE3& inVec2) {
        const __m128 testResult = _mm_cmpgt_ps(inVec1.vec, inVec2.vec);
        const __m128 ones       = _mm_set1_ps(1);
        return _mm_and_ps(testResult, ones);
    }

    inline static InaVecSSE3 IsEqual(const InaVecSSE3& inVec1, const InaVecSSE3& inVec2) {
        const __m128 testResult = _mm_cmpeq_ps(inVec1.vec, inVec2.vec);
        const __m128 ones       = _mm_set1_ps(1);
        return _mm_and_ps(testResult, ones);
    }

    inline static InaVecSSE3 IsNotEqual(const InaVecSSE3& inVec1, const InaVecSSE3& inVec2) {
        const __m128 testResult = _mm_cmpneq_ps(inVec1.vec, inVec2.vec);
        const __m128 ones       = _mm_set1_ps(1);
        return _mm_and_ps(testResult, ones);
    }

    inline static InaVecMaskSSE3< float > IsLowerOrEqualMask(const InaVecSSE3& inVec1, const InaVecSSE3& inVec2) {
        return _mm_castps_si128(_mm_cmple_ps(inVec1.vec, inVec2.vec));
    }

    inline static InaVecMaskSSE3< float > IsLowerMask(const InaVecSSE3& inVec1, const InaVecSSE3& inVec2) {
        return _mm_castps_si128(_mm_cmplt_ps(inVec1.vec, inVec2.vec));
    }

    inline static InaVecMaskSSE3< float > IsGreaterOrEqualMask(const InaVecSSE3& inVec1, const InaVecSSE3& inVec2) {
        return _mm_castps_si128(_mm_cmpge_ps(inVec1.vec, inVec2.vec));
    }

    inline static InaVecMaskSSE3< float > IsGreaterMask(const InaVecSSE3& inVec1, const InaVecSSE3& inVec2) {
        return _mm_castps_si128(_mm_cmpgt_ps(inVec1.vec, inVec2.vec));
    }

    inline static InaVecMaskSSE3< float > IsEqualMask(const InaVecSSE3& inVec1, const InaVecSSE3& inVec2) {
        return _mm_castps_si128(_mm_cmpeq_ps(inVec1.vec, inVec2.vec));
    }

    inline static InaVecMaskSSE3< float > IsNotEqualMask(const InaVecSSE3& inVec1, const InaVecSSE3& inVec2) {
        return _mm_castps_si128(_mm_cmpneq_ps(inVec1.vec, inVec2.vec));
    }

    inline static InaVecSSE3 BitsAnd(const InaVecSSE3& inVec1, const InaVecSSE3& inVec2) {
        return _mm_and_ps(inVec1.vec, inVec2.vec);
    }

    inline static InaVecSSE3 BitsNotAnd(const InaVecSSE3& inVec1, const InaVecSSE3& inVec2) {
        return _mm_andnot_ps(inVec1.vec, inVec2.vec);
    }

    inline static InaVecSSE3 BitsOr(const InaVecSSE3& inVec1, const InaVecSSE3& inVec2) {
        return _mm_or_ps(inVec1.vec, inVec2.vec);
    }

    inline static InaVecSSE3 BitsXor(const InaVecSSE3& inVec1, const InaVecSSE3& inVec2) {
        return _mm_xor_ps(inVec1.vec, inVec2.vec);
    }

    inline static const char* GetName() {
        return "InaVecSSE3<float>";
    }

    inline static InaIfElse< InaVecSSE3< float > >::ThenClass If(const InaVecMaskSSE3< float >& inTest) {
        return InaIfElse< InaVecSSE3< float > >::IfClass().If(inTest);
    }

    inline static InaVecSSE3 IfElse(const InaVecMaskSSE3< float >& inMask, const InaVecSSE3& inIfTrue, const InaVecSSE3& inIfFalse) {
        return _mm_or_ps(IfTrue(inMask, inIfTrue.vec).vec,
                         IfFalse(inMask, inIfFalse.vec).vec);
    }

    inline static InaVecSSE3 IfTrue(const InaVecMaskSSE3< float >& inMask, const InaVecSSE3& inIfTrue) {
        return _mm_and_ps(_mm_castsi128_ps(inMask.getMask()), inIfTrue.vec);
    }

    inline static InaVecSSE3 IfFalse(const InaVecMaskSSE3< float >& inMask, const InaVecSSE3& inIfFalse) {
        return _mm_andnot_ps(_mm_castsi128_ps(inMask.getMask()), inIfFalse.vec);
    }

    // Inner operators
    inline InaVecSSE3< float >& operator+=(const InaVecSSE3< float >& inVec) {
        vec = _mm_add_ps(vec, inVec.vec);
        return *this;
    }

    inline InaVecSSE3< float >& operator-=(const InaVecSSE3< float >& inVec) {
        vec = _mm_sub_ps(vec, inVec.vec);
        return *this;
    }

    inline InaVecSSE3< float >& operator/=(const InaVecSSE3< float >& inVec) {
        vec = _mm_div_ps(vec, inVec.vec);
        return *this;
    }

    inline InaVecSSE3< float >& operator*=(const InaVecSSE3< float >& inVec) {
        vec = _mm_mul_ps(vec, inVec.vec);
        return *this;
    }

    inline InaVecSSE3< float > operator-() const {
        const __m128 minus0 = _mm_castsi128_ps(_mm_set1_epi32(static_cast< int >(0x80000000)));
        return _mm_xor_ps(vec, minus0);
    }

    inline InaVecSSE3< float > pow(std::size_t power) const {
        return InaUtils::FastPow< InaVecSSE3< float > >(*this, power);
    }
};

// Bits operators
inline InaVecSSE3< float > operator&(const InaVecSSE3< float >& inVec1, const InaVecSSE3< float >& inVec2) {
    return InaVecSSE3< float >::BitsAnd(inVec1, inVec2);
}

inline InaVecSSE3< float > operator|(const InaVecSSE3< float >& inVec1, const InaVecSSE3< float >& inVec2) {
    return InaVecSSE3< float >::BitsOr(inVec1, inVec2);
}

inline InaVecSSE3< float > operator^(const InaVecSSE3< float >& inVec1, const InaVecSSE3< float >& inVec2) {
    return InaVecSSE3< float >::BitsXor(inVec1, inVec2);
}

// Dual operators
inline InaVecSSE3< float > operator+(const InaVecSSE3< float >& inVec1, const InaVecSSE3< float >& inVec2) {
    return _mm_add_ps(inVec1.getVec(), inVec2.getVec());
}

inline InaVecSSE3< float > operator-(const InaVecSSE3< float >& inVec1, const InaVecSSE3< float >& inVec2) {
    return _mm_sub_ps(inVec1.getVec(), inVec2.getVec());
}

inline InaVecSSE3< float > operator/(const InaVecSSE3< float >& inVec1, const InaVecSSE3< float >& inVec2) {
    return _mm_div_ps(inVec1.getVec(), inVec2.getVec());
}

inline InaVecSSE3< float > operator*(const InaVecSSE3< float >& inVec1, const InaVecSSE3< float >& inVec2) {
    return _mm_mul_ps(inVec1.getVec(), inVec2.getVec());
}

// Tests and comparions
inline InaVecMaskSSE3< float > operator<(const InaVecSSE3< float >& inVec1, const InaVecSSE3< float >& inVec2) {
    return InaVecSSE3< float >::IsLowerMask(inVec1, inVec2);
}

inline InaVecMaskSSE3< float > operator<=(const InaVecSSE3< float >& inVec1, const InaVecSSE3< float >& inVec2) {
    return InaVecSSE3< float >::IsLowerOrEqualMask(inVec1, inVec2);
}

inline InaVecMaskSSE3< float > operator>(const InaVecSSE3< float >& inVec1, const InaVecSSE3< float >& inVec2) {
    return InaVecSSE3< float >::IsGreaterMask(inVec1, inVec2);
}

inline InaVecMaskSSE3< float > operator>=(const InaVecSSE3< float >& inVec1, const InaVecSSE3< float >& inVec2) {
    return InaVecSSE3< float >::IsGreaterOrEqualMask(inVec1, inVec2);
}

inline InaVecMaskSSE3< float > operator==(const InaVecSSE3< float >& inVec1, const InaVecSSE3< float >& inVec2) {
    return InaVecSSE3< float >::IsEqualMask(inVec1, inVec2);
}

inline InaVecMaskSSE3< float > operator!=(const InaVecSSE3< float >& inVec1, const InaVecSSE3< float >& inVec2) {
    return InaVecSSE3< float >::IsNotEqualMask(inVec1, inVec2);
}


#endif
