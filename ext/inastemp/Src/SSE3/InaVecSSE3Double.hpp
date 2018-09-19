///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef INAVECSSE3DOUBLE_HPP
#define INAVECSSE3DOUBLE_HPP

#include "InastempConfig.h"
#include "Common/InaIfElse.hpp"
#include "Common/InaUtils.hpp"

#ifndef INASTEMP_USE_SSE3
#error InaVecSSE3<double> is included but SSE3 is not enable in the configuration
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
class alignas(16) InaVecMaskSSE3< double > {
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

    explicit operator __m128i() const {
        return mask;
    }

    __m128i getMask() const {
        return mask;
    }

    // Bool data type compatibility
    inline explicit InaVecMaskSSE3(const bool inBool) {
        mask = (inBool ? _mm_set1_epi64x(static_cast< long long >(0xFFFFFFFFFFFFFFFFL)) : _mm_setzero_si128());
    }

    inline InaVecMaskSSE3& operator=(const bool inBool) {
        mask = (inBool ? _mm_set1_epi64x(static_cast< long long >(0xFFFFFFFFFFFFFFFFL)) : _mm_setzero_si128());
        return (*this);
    }

    // Binary methods
    inline InaVecMaskSSE3 Not() const {
        return NotAnd(mask, _mm_set1_epi64x(static_cast< long long >(0xFFFFFFFFFFFFFFFFL)));
    }

    inline bool isAllTrue() const {
        // true if all FF
        return _mm_movemask_epi8(_mm_cmpeq_epi32(mask, _mm_set1_epi64x(static_cast< long long >(0xFFFFFFFFFFFFFFFFL)))) == 0xFFFF;
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
inline InaVecMaskSSE3< double > operator&(const InaVecMaskSSE3< double >& inMask1, const InaVecMaskSSE3< double >& inMask2) {
    return InaVecMaskSSE3< double >::And(inMask1, inMask2);
}

inline InaVecMaskSSE3< double > operator|(const InaVecMaskSSE3< double >& inMask1, const InaVecMaskSSE3< double >& inMask2) {
    return InaVecMaskSSE3< double >::Or(inMask1, inMask2);
}

inline InaVecMaskSSE3< double > operator^(const InaVecMaskSSE3< double >& inMask1, const InaVecMaskSSE3< double >& inMask2) {
    return InaVecMaskSSE3< double >::Xor(inMask1, inMask2);
}

inline bool operator==(const InaVecMaskSSE3< double >& inMask1, const InaVecMaskSSE3< double >& inMask2) {
    return InaVecMaskSSE3< double >::IsEqual(inMask1, inMask2);
}

inline bool operator!=(const InaVecMaskSSE3< double >& inMask1, const InaVecMaskSSE3< double >& inMask2) {
    return InaVecMaskSSE3< double >::IsNotEqual(inMask1, inMask2);
}

// Vec type
template <>
class alignas(16) InaVecSSE3< double > {
protected:
    __m128d vec;

public:
    using VecRawType            = __m128d;
    using MaskType              = InaVecMaskSSE3< double >;
    using RealType              = double;
    static const int VecLength  = 2;
    static const int Alignement = 16;

    inline InaVecSSE3() {
    }
    inline InaVecSSE3(const InaVecSSE3&) = default;
    inline InaVecSSE3& operator=(const InaVecSSE3&) = default;

    // Constructor from raw type
    inline /*not explicit*/ InaVecSSE3(const __m128d inVec)
    : vec(inVec) {
    }

    inline InaVecSSE3& operator=(const __m128d inVec) {
        vec = inVec;
        return *this;
    }

    inline void setFromRawType(const __m128d inVec) {
        vec = inVec;
    }

    inline explicit operator __m128d() const {
        return vec;
    }

    inline __m128d getVec() const {
        return vec;
    }

    // Constructor from scalar
    inline /*not explicit*/ InaVecSSE3(const double inVal)
    : vec(_mm_set1_pd(inVal)) {
    }

    inline InaVecSSE3& operator=(const double inVal) {
        vec = _mm_set1_pd(inVal);
        return *this;
    }

    inline void setFromScalar(const double inVal) {
        vec = _mm_set1_pd(inVal);
    }

    // Constructor from vec
    inline InaVecSSE3(const std::initializer_list< double > lst)
    : InaVecSSE3(lst.begin()) {
    }

    inline explicit InaVecSSE3(const double ptr[])
    : vec(_mm_loadu_pd(ptr)) {
    }

    inline InaVecSSE3& setFromArray(const double ptr[]) {
        vec = _mm_loadu_pd(ptr);
        return *this;
    }

    inline InaVecSSE3& setFromAlignedArray(const double ptr[]) {
        vec = _mm_load_pd(ptr);
        return *this;
    }

    inline InaVecSSE3& setFromIndirectArray(const double values[], const int inIndirection[]) {
        vec = _mm_set_pd(
            values[inIndirection[1]],
            values[inIndirection[0]]);
        return *this;
    }

    inline InaVecSSE3& setFromIndirect2DArray(const double inArray[], const int inIndirection1[],
                                              const int inLeadingDimension, const int inIndirection2[]) {
        vec = _mm_set_pd(
            inArray[inIndirection1[1] * inLeadingDimension + inIndirection2[1]],
            inArray[inIndirection1[0] * inLeadingDimension + inIndirection2[0]]);
        return *this;
    }

    // Move back to array
    inline void storeInArray(double ptr[]) const {
        _mm_storeu_pd(ptr, vec);
    }

    inline void storeInAlignedArray(double ptr[]) const {
        _mm_store_pd(ptr, vec);
    }

    // Acce to individual values
    inline double at(const int index) const {
        alignas(Alignement) double allval[VecLength];
        _mm_store_pd(allval, vec);
        return allval[index];
    }

    // Horizontal operation
    inline double horizontalSum() const {
        const __m128d res = _mm_add_pd(vec, _mm_shuffle_pd(vec, vec, 1));
        return _mm_cvtsd_f64(res);
    }

    inline double horizontalMul() const {
        const __m128d res = _mm_mul_pd(vec, _mm_shuffle_pd(vec, vec, 1));
        return _mm_cvtsd_f64(res);
    }

    inline InaVecSSE3 sqrt() const {
        return _mm_sqrt_pd(vec);
    }

    inline InaVecSSE3 exp() const {
#ifdef __INTEL_COMPILER
        return _mm_exp_pd(vec);
#else
        const __m128d COEFF_LOG2E = _mm_set1_pd(double(InaFastExp::CoeffLog2E()));
        const __m128d COEFF_A     = _mm_set1_pd(double(InaFastExp::CoeffA64()));
        const __m128d COEFF_B     = _mm_set1_pd(double(InaFastExp::CoeffB64()));
        const __m128d COEFF_P5_X  = _mm_set1_pd(double(InaFastExp::GetCoefficient9_8()));
        const __m128d COEFF_P5_Y  = _mm_set1_pd(double(InaFastExp::GetCoefficient9_7()));
        const __m128d COEFF_P5_Z  = _mm_set1_pd(double(InaFastExp::GetCoefficient9_6()));
        const __m128d COEFF_P5_A  = _mm_set1_pd(double(InaFastExp::GetCoefficient9_5()));
        const __m128d COEFF_P5_B  = _mm_set1_pd(double(InaFastExp::GetCoefficient9_4()));
        const __m128d COEFF_P5_C  = _mm_set1_pd(double(InaFastExp::GetCoefficient9_3()));
        const __m128d COEFF_P5_D  = _mm_set1_pd(double(InaFastExp::GetCoefficient9_2()));
        const __m128d COEFF_P5_E  = _mm_set1_pd(double(InaFastExp::GetCoefficient9_1()));
        const __m128d COEFF_P5_F  = _mm_set1_pd(double(InaFastExp::GetCoefficient9_0()));

        __m128d x = _mm_mul_pd(vec, COEFF_LOG2E);

        const __m128d fractional_part = _mm_sub_pd(x, InaVecSSE3(x).floor().vec);

        __m128d factor = _mm_add_pd(_mm_mul_pd(_mm_add_pd(
                                                   _mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(
                                                                                        _mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(
                                                                                                                             _mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(
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

        x = _mm_sub_pd(x, factor);

        x = _mm_add_pd(_mm_mul_pd(COEFF_A, x), COEFF_B);

        alignas(64) long int allvalint[VecLength] = { _mm_cvtsd_si64(x),
                                                      _mm_cvtsd_si64(_mm_shuffle_pd(x, x, 1)) };

        return _mm_castsi128_pd(_mm_set_epi64x(allvalint[1], allvalint[0]));
#endif
    }

    inline InaVecSSE3 expLowAcc() const {
        const __m128d COEFF_LOG2E = _mm_set1_pd(double(InaFastExp::CoeffLog2E()));
        const __m128d COEFF_A     = _mm_set1_pd(double(InaFastExp::CoeffA64()));
        const __m128d COEFF_B     = _mm_set1_pd(double(InaFastExp::CoeffB64()));
        const __m128d COEFF_P5_C  = _mm_set1_pd(double(InaFastExp::GetCoefficient4_3()));
        const __m128d COEFF_P5_D  = _mm_set1_pd(double(InaFastExp::GetCoefficient4_2()));
        const __m128d COEFF_P5_E  = _mm_set1_pd(double(InaFastExp::GetCoefficient4_1()));
        const __m128d COEFF_P5_F  = _mm_set1_pd(double(InaFastExp::GetCoefficient4_0()));

        __m128d x = _mm_mul_pd(vec, COEFF_LOG2E);

        const __m128d fractional_part = _mm_sub_pd(x, InaVecSSE3(x).floor().vec);

        __m128d factor = _mm_add_pd(_mm_mul_pd(_mm_add_pd(
                                                   _mm_mul_pd(_mm_add_pd(_mm_mul_pd(
                                                                             COEFF_P5_C, fractional_part),
                                                                         COEFF_P5_D),
                                                              fractional_part),
                                                   COEFF_P5_E),
                                               fractional_part),
                                    COEFF_P5_F);

        x = _mm_sub_pd(x, factor);

        x = _mm_add_pd(_mm_mul_pd(COEFF_A, x), COEFF_B);

        alignas(64) long int allvalint[VecLength] = { _mm_cvtsd_si64(x),
                                                      _mm_cvtsd_si64(_mm_shuffle_pd(x, x, 1)) };

        return _mm_castsi128_pd(_mm_set_epi64x(allvalint[1], allvalint[0]));
    }

    inline InaVecSSE3 rsqrt() const {
        return _mm_set1_pd(1) / _mm_sqrt_pd(vec);
    }

    inline InaVecSSE3 abs() const {
        const __m128d minus0 = _mm_castsi128_pd(_mm_set1_epi64x(static_cast< long long >(0x8000000000000000L)));
        return _mm_andnot_pd(minus0, vec);
    }

    inline InaVecSSE3 floor() const {
        alignas(Alignement) double allval[VecLength];
        _mm_store_pd(allval, vec);
        for(int idx = 0; idx < VecLength; ++idx) {
            allval[idx] = std::floor(allval[idx]);
        }
        return _mm_loadu_pd(allval);
    }

    inline InaVecSSE3 signOf() const {
        const __m128d minus0 = _mm_castsi128_pd(_mm_set1_epi64x(static_cast< long long >(0x8000000000000000L)));
        const __m128d signs  = _mm_and_pd(vec, minus0);
        return _mm_andnot_pd(_mm_cmpeq_pd(_mm_setzero_pd(), vec), _mm_or_pd(signs, _mm_set1_pd(1)));
    }

    inline InaVecSSE3 isPositive() const {
        const __m128d greater = _mm_cmple_pd(_mm_setzero_pd(), vec);
        const __m128d ones    = _mm_set1_pd(1);
        return _mm_and_pd(greater, ones);
    }

    inline InaVecSSE3 isNegative() const {
        const __m128d less = _mm_cmpge_pd(_mm_setzero_pd(), vec);
        const __m128d ones = _mm_set1_pd(1);
        return _mm_and_pd(less, ones);
    }

    inline InaVecSSE3 isPositiveStrict() const {
        const __m128d greater = _mm_cmplt_pd(_mm_setzero_pd(), vec);
        const __m128d ones    = _mm_set1_pd(1);
        return _mm_and_pd(greater, ones);
    }

    inline InaVecSSE3 isNegativeStrict() const {
        const __m128d less = _mm_cmpgt_pd(_mm_setzero_pd(), vec);
        const __m128d ones = _mm_set1_pd(1);
        return _mm_and_pd(less, ones);
    }

    inline InaVecSSE3 isZero() const {
        const __m128d equalZero = _mm_cmpeq_pd(_mm_setzero_pd(), vec);
        const __m128d ones      = _mm_set1_pd(1);
        return _mm_and_pd(equalZero, ones);
    }

    inline InaVecSSE3 isNotZero() const {
        const __m128d equalZero = _mm_cmpeq_pd(_mm_setzero_pd(), vec);
        const __m128d ones      = _mm_set1_pd(1);
        return _mm_andnot_pd(equalZero, ones);
    }

    inline InaVecMaskSSE3< double > isPositiveMask() const {
        return _mm_castpd_si128(_mm_cmple_pd(_mm_setzero_pd(), vec));
    }

    inline InaVecMaskSSE3< double > isNegativeMask() const {
        return _mm_castpd_si128(_mm_cmpge_pd(_mm_setzero_pd(), vec));
    }

    inline InaVecMaskSSE3< double > isPositiveStrictMask() const {
        return _mm_castpd_si128(_mm_cmplt_pd(_mm_setzero_pd(), vec));
    }

    inline InaVecMaskSSE3< double > isNegativeStrictMask() const {
        return _mm_castpd_si128(_mm_cmpgt_pd(_mm_setzero_pd(), vec));
    }

    inline InaVecMaskSSE3< double > isZeroMask() const {
        return _mm_castpd_si128(_mm_cmpeq_pd(_mm_setzero_pd(), vec));
    }

    inline InaVecMaskSSE3< double > isNotZeroMask() const {
        return _mm_castpd_si128(_mm_cmpneq_pd(_mm_setzero_pd(), vec));
    }

    // Static basic methods
    inline static InaVecSSE3 GetZero() {
        return InaVecSSE3(_mm_setzero_pd());
    }

    inline static InaVecSSE3 GetOne() {
        return InaVecSSE3(_mm_set1_pd(1));
    }

    inline static InaVecSSE3 Min(const InaVecSSE3& inVec1, const InaVecSSE3& inVec2) {
        return _mm_min_pd(inVec1.vec, inVec2.vec);
    }

    inline static InaVecSSE3 Max(const InaVecSSE3& inVec1, const InaVecSSE3& inVec2) {
        return _mm_max_pd(inVec1.vec, inVec2.vec);
    }

    inline static InaVecSSE3 IsLowerOrEqual(const InaVecSSE3& inVec1, const InaVecSSE3& inVec2) {
        const __m128d testResult = _mm_cmple_pd(inVec1.vec, inVec2.vec);
        const __m128d ones       = _mm_set1_pd(1);
        return _mm_and_pd(testResult, ones);
    }

    inline static InaVecSSE3 IsLower(const InaVecSSE3& inVec1, const InaVecSSE3& inVec2) {
        const __m128d testResult = _mm_cmplt_pd(inVec1.vec, inVec2.vec);
        const __m128d ones       = _mm_set1_pd(1);
        return _mm_and_pd(testResult, ones);
    }

    inline static InaVecSSE3 IsGreaterOrEqual(const InaVecSSE3& inVec1, const InaVecSSE3& inVec2) {
        const __m128d testResult = _mm_cmpge_pd(inVec1.vec, inVec2.vec);
        const __m128d ones       = _mm_set1_pd(1);
        return _mm_and_pd(testResult, ones);
    }

    inline static InaVecSSE3 IsGreater(const InaVecSSE3& inVec1, const InaVecSSE3& inVec2) {
        const __m128d testResult = _mm_cmpgt_pd(inVec1.vec, inVec2.vec);
        const __m128d ones       = _mm_set1_pd(1);
        return _mm_and_pd(testResult, ones);
    }

    inline static InaVecSSE3 IsEqual(const InaVecSSE3& inVec1, const InaVecSSE3& inVec2) {
        const __m128d testResult = _mm_cmpeq_pd(inVec1.vec, inVec2.vec);
        const __m128d ones       = _mm_set1_pd(1);
        return _mm_and_pd(testResult, ones);
    }

    inline static InaVecSSE3 IsNotEqual(const InaVecSSE3& inVec1, const InaVecSSE3& inVec2) {
        const __m128d testResult = _mm_cmpneq_pd(inVec1.vec, inVec2.vec);
        const __m128d ones       = _mm_set1_pd(1);
        return _mm_and_pd(testResult, ones);
    }

    inline static InaVecMaskSSE3< double > IsLowerOrEqualMask(const InaVecSSE3& inVec1, const InaVecSSE3& inVec2) {
        return _mm_castpd_si128(_mm_cmple_pd(inVec1.vec, inVec2.vec));
    }

    inline static InaVecMaskSSE3< double > IsLowerMask(const InaVecSSE3& inVec1, const InaVecSSE3& inVec2) {
        return _mm_castpd_si128(_mm_cmplt_pd(inVec1.vec, inVec2.vec));
    }

    inline static InaVecMaskSSE3< double > IsGreaterOrEqualMask(const InaVecSSE3& inVec1, const InaVecSSE3& inVec2) {
        return _mm_castpd_si128(_mm_cmpge_pd(inVec1.vec, inVec2.vec));
    }

    inline static InaVecMaskSSE3< double > IsGreaterMask(const InaVecSSE3& inVec1, const InaVecSSE3& inVec2) {
        return _mm_castpd_si128(_mm_cmpgt_pd(inVec1.vec, inVec2.vec));
    }

    inline static InaVecMaskSSE3< double > IsEqualMask(const InaVecSSE3& inVec1, const InaVecSSE3& inVec2) {
        return _mm_castpd_si128(_mm_cmpeq_pd(inVec1.vec, inVec2.vec));
    }

    inline static InaVecMaskSSE3< double > IsNotEqualMask(const InaVecSSE3& inVec1, const InaVecSSE3& inVec2) {
        return _mm_castpd_si128(_mm_cmpneq_pd(inVec1.vec, inVec2.vec));
    }

    inline static InaVecSSE3 BitsAnd(const InaVecSSE3& inVec1, const InaVecSSE3& inVec2) {
        return _mm_and_pd(inVec1.vec, inVec2.vec);
    }

    inline static InaVecSSE3 BitsNotAnd(const InaVecSSE3& inVec1, const InaVecSSE3& inVec2) {
        return _mm_andnot_pd(inVec1.vec, inVec2.vec);
    }

    inline static InaVecSSE3 BitsOr(const InaVecSSE3& inVec1, const InaVecSSE3& inVec2) {
        return _mm_or_pd(inVec1.vec, inVec2.vec);
    }

    inline static InaVecSSE3 BitsXor(const InaVecSSE3& inVec1, const InaVecSSE3& inVec2) {
        return _mm_xor_pd(inVec1.vec, inVec2.vec);
    }

    inline static const char* GetName() {
        return "InaVecSSE3<double>";
    }

    inline static InaIfElse< InaVecSSE3< double > >::ThenClass If(const InaVecMaskSSE3< double >& inTest) {
        return InaIfElse< InaVecSSE3< double > >::IfClass().If(inTest);
    }

    inline static InaVecSSE3 IfElse(const InaVecMaskSSE3< double >& inMask, const InaVecSSE3& inIfTrue, const InaVecSSE3& inIfFalse) {
        return _mm_or_pd(IfTrue(inMask, inIfTrue.vec).vec,
                         IfFalse(inMask, inIfFalse.vec).vec);
    }

    inline static InaVecSSE3 IfTrue(const InaVecMaskSSE3< double >& inMask, const InaVecSSE3& inIfTrue) {
        return _mm_and_pd(_mm_castsi128_pd(inMask.getMask()), inIfTrue.vec);
    }

    inline static InaVecSSE3 IfFalse(const InaVecMaskSSE3< double >& inMask, const InaVecSSE3& inIfFalse) {
        return _mm_andnot_pd(_mm_castsi128_pd(inMask.getMask()), inIfFalse.vec);
    }

    // Inner operators
    inline InaVecSSE3< double >& operator+=(const InaVecSSE3< double >& inVec) {
        vec = _mm_add_pd(vec, inVec.vec);
        return *this;
    }

    inline InaVecSSE3< double >& operator-=(const InaVecSSE3< double >& inVec) {
        vec = _mm_sub_pd(vec, inVec.vec);
        return *this;
    }

    inline InaVecSSE3< double >& operator/=(const InaVecSSE3< double >& inVec) {
        vec = _mm_div_pd(vec, inVec.vec);
        return *this;
    }

    inline InaVecSSE3< double >& operator*=(const InaVecSSE3< double >& inVec) {
        vec = _mm_mul_pd(vec, inVec.vec);
        return *this;
    }

    inline InaVecSSE3< double > operator-() const {
        const __m128d minus0 = _mm_castsi128_pd(_mm_set1_epi64x(static_cast< long long >(0x8000000000000000L)));
        return _mm_xor_pd(vec, minus0);
    }

    inline InaVecSSE3< double > pow(std::size_t power) const {
        return InaUtils::FastPow< InaVecSSE3< double > >(*this, power);
    }
};

// Bits operators
inline InaVecSSE3< double > operator&(const InaVecSSE3< double >& inVec1, const InaVecSSE3< double >& inVec2) {
    return InaVecSSE3< double >::BitsAnd(inVec1, inVec2);
}

inline InaVecSSE3< double > operator|(const InaVecSSE3< double >& inVec1, const InaVecSSE3< double >& inVec2) {
    return InaVecSSE3< double >::BitsOr(inVec1, inVec2);
}

inline InaVecSSE3< double > operator^(const InaVecSSE3< double >& inVec1, const InaVecSSE3< double >& inVec2) {
    return InaVecSSE3< double >::BitsXor(inVec1, inVec2);
}

// Dual operators
inline InaVecSSE3< double > operator+(const InaVecSSE3< double >& inVec1, const InaVecSSE3< double >& inVec2) {
    return _mm_add_pd(inVec1.getVec(), inVec2.getVec());
}

inline InaVecSSE3< double > operator-(const InaVecSSE3< double >& inVec1, const InaVecSSE3< double >& inVec2) {
    return _mm_sub_pd(inVec1.getVec(), inVec2.getVec());
}

inline InaVecSSE3< double > operator/(const InaVecSSE3< double >& inVec1, const InaVecSSE3< double >& inVec2) {
    return _mm_div_pd(inVec1.getVec(), inVec2.getVec());
}

inline InaVecSSE3< double > operator*(const InaVecSSE3< double >& inVec1, const InaVecSSE3< double >& inVec2) {
    return _mm_mul_pd(inVec1.getVec(), inVec2.getVec());
}

// Tests and comparions
inline InaVecMaskSSE3< double > operator<(const InaVecSSE3< double >& inVec1, const InaVecSSE3< double >& inVec2) {
    return InaVecSSE3< double >::IsLowerMask(inVec1, inVec2);
}

inline InaVecMaskSSE3< double > operator<=(const InaVecSSE3< double >& inVec1, const InaVecSSE3< double >& inVec2) {
    return InaVecSSE3< double >::IsLowerOrEqualMask(inVec1, inVec2);
}

inline InaVecMaskSSE3< double > operator>(const InaVecSSE3< double >& inVec1, const InaVecSSE3< double >& inVec2) {
    return InaVecSSE3< double >::IsGreaterMask(inVec1, inVec2);
}

inline InaVecMaskSSE3< double > operator>=(const InaVecSSE3< double >& inVec1, const InaVecSSE3< double >& inVec2) {
    return InaVecSSE3< double >::IsGreaterOrEqualMask(inVec1, inVec2);
}

inline InaVecMaskSSE3< double > operator==(const InaVecSSE3< double >& inVec1, const InaVecSSE3< double >& inVec2) {
    return InaVecSSE3< double >::IsEqualMask(inVec1, inVec2);
}

inline InaVecMaskSSE3< double > operator!=(const InaVecSSE3< double >& inVec1, const InaVecSSE3< double >& inVec2) {
    return InaVecSSE3< double >::IsNotEqualMask(inVec1, inVec2);
}

#endif
