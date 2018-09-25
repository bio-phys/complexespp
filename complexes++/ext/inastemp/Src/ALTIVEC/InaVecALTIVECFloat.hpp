///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef INAVECALTIVECFLOAT_HPP
#define INAVECALTIVECFLOAT_HPP

#include "InastempConfig.h"
#include "Common/InaIfElse.hpp"
#include "Common/InaUtils.hpp"

#ifndef INASTEMP_USE_ALTIVEC
#error InaVecALTIVEC<float> is included but ALTIVEC is not enable in the configuration
#endif

#include "Common/InaFastExp.hpp"

#include <altivec.h>
#undef bool
#undef vector
#undef pixel

#include <initializer_list>

// Forward declarations
template < class RealType >
class InaVecMaskALTIVEC;

template < class RealType >
class InaVecALTIVEC;


// Mask type
template <>
class alignas(16) InaVecMaskALTIVEC< float > {
    __vector __bool int mask;

public:
    // Classic constructors
    inline InaVecMaskALTIVEC() {
    }

    InaVecMaskALTIVEC(const InaVecMaskALTIVEC&) = default;
    inline InaVecMaskALTIVEC& operator=(const InaVecMaskALTIVEC&) = default;

    // Native data type compatibility
    inline /*not explicit*/ InaVecMaskALTIVEC(const __vector __bool int inMask)
    : mask(inMask) {
    }

    inline InaVecMaskALTIVEC& operator=(const __vector __bool int inMask) {
        mask = inMask;
        return (*this);
    }

    inline explicit operator __vector __bool int() const {
        return mask;
    }

    inline __vector __bool int getMask() const {
        return mask;
    }

    // Bool data type compatibility
    inline explicit InaVecMaskALTIVEC(const bool inBool) {
        const __vector __bool int tmpMaskFF = reinterpret_cast< __vector __bool int >(vec_splats(0xFFFFFFFFU));
        mask                                = (inBool ? tmpMaskFF : reinterpret_cast< __vector __bool int >(vec_splats(0x0U)));
    }

    inline InaVecMaskALTIVEC& operator=(const bool inBool) {
        const __vector __bool int tmpMaskFF = reinterpret_cast< __vector __bool int >(vec_splats(0xFFFFFFFFU));
        mask                                = (inBool ? tmpMaskFF : reinterpret_cast< __vector __bool int >(vec_splats(0x0U)));
        return (*this);
    }

    // Binary methods
    inline InaVecMaskALTIVEC Not() const {
        const __vector __bool int tmpMaskFF = reinterpret_cast< __vector __bool int >(vec_splats(0xFFFFFFFFU));
        return NotAnd(mask, tmpMaskFF);
    }

    inline bool isAllTrue() const {
        // true if all FF => !FF => 0 & FF => 0
        const __vector __bool int tmpMaskFF = reinterpret_cast< __vector __bool int >(vec_splats(0xFFFFFFFFU));
        ;
        const int res = vec_all_eq(mask, tmpMaskFF);
        return static_cast< bool >(res);
    }

    inline bool isAllFalse() const {
        // true if all zero
        const int res = vec_all_eq(mask, vec_xor(mask, mask));
        return static_cast< bool >(res);
    }

    // Double args methods
    inline static InaVecMaskALTIVEC And(const InaVecMaskALTIVEC& inMask1, const InaVecMaskALTIVEC& inMask2) {
        return InaVecMaskALTIVEC(vec_and(inMask1.mask, inMask2.mask));
    }

    inline static InaVecMaskALTIVEC NotAnd(const InaVecMaskALTIVEC& inMask1, const InaVecMaskALTIVEC& inMask2) {
        return InaVecMaskALTIVEC(reinterpret_cast< __vector __bool int >(
            vec_nand(reinterpret_cast< __vector unsigned int >(inMask1.mask),
                     reinterpret_cast< __vector unsigned int >(inMask2.mask))));
    }

    inline static InaVecMaskALTIVEC Or(const InaVecMaskALTIVEC& inMask1, const InaVecMaskALTIVEC& inMask2) {
        return InaVecMaskALTIVEC(vec_or(inMask1.mask, inMask2.mask));
    }

    inline static InaVecMaskALTIVEC Xor(const InaVecMaskALTIVEC& inMask1, const InaVecMaskALTIVEC& inMask2) {
        return InaVecMaskALTIVEC(vec_xor(inMask1.mask, inMask2.mask));
    }

    inline static bool IsEqual(const InaVecMaskALTIVEC& inMask1, const InaVecMaskALTIVEC& inMask2) {
        const int res = vec_all_eq(inMask1.mask, inMask2.mask);
        return static_cast< bool >(res);
    }

    inline static bool IsNotEqual(const InaVecMaskALTIVEC& inMask1, const InaVecMaskALTIVEC& inMask2) {
        const int res = !vec_all_eq(inMask1.mask, inMask2.mask);
        return static_cast< bool >(res);
    }
};

// Mask must have operators
inline InaVecMaskALTIVEC< float > operator&(const InaVecMaskALTIVEC< float >& inMask1, const InaVecMaskALTIVEC< float >& inMask2) {
    return InaVecMaskALTIVEC< float >::And(inMask1, inMask2);
}

inline InaVecMaskALTIVEC< float > operator|(const InaVecMaskALTIVEC< float >& inMask1, const InaVecMaskALTIVEC< float >& inMask2) {
    return InaVecMaskALTIVEC< float >::Or(inMask1, inMask2);
}

inline InaVecMaskALTIVEC< float > operator^(const InaVecMaskALTIVEC< float >& inMask1, const InaVecMaskALTIVEC< float >& inMask2) {
    return InaVecMaskALTIVEC< float >::Xor(inMask1, inMask2);
}

inline bool operator==(const InaVecMaskALTIVEC< float >& inMask1, const InaVecMaskALTIVEC< float >& inMask2) {
    return InaVecMaskALTIVEC< float >::IsEqual(inMask1, inMask2);
}

inline bool operator!=(const InaVecMaskALTIVEC< float >& inMask1, const InaVecMaskALTIVEC< float >& inMask2) {
    return InaVecMaskALTIVEC< float >::IsNotEqual(inMask1, inMask2);
}

// Vec type
template <>
class alignas(16) InaVecALTIVEC< float > {
protected:
    __vector float vec;

public:
    using VecRawType            = __vector float;
    using MaskType              = InaVecMaskALTIVEC< float >;
    using RealType              = float;
    static const int VecLength  = 4;
    static const int Alignement = 16;

    inline InaVecALTIVEC() {
    }
    inline InaVecALTIVEC(const InaVecALTIVEC&) = default;
    inline InaVecALTIVEC& operator=(const InaVecALTIVEC&) = default;

    // Constructor from raw type
    inline /*not explicit*/ InaVecALTIVEC(const __vector float inVec)
    : vec(inVec) {
    }

    inline InaVecALTIVEC& operator=(const __vector float inVec) {
        vec = inVec;
        return *this;
    }

    inline void setFromRawType(const __vector float inVec) {
        vec = inVec;
    }

    inline explicit operator __vector float() const {
        return vec;
    }

    inline __vector float getVec() const {
        return vec;
    }

    // Constructor from scalar
    inline /*not explicit*/ InaVecALTIVEC(const float val)
    : vec(vec_splats(val)) {
    }

    inline InaVecALTIVEC& operator=(const float val) {
        vec = vec_splats(val);
        return *this;
    }

    // Constructor from vec
    inline InaVecALTIVEC(const std::initializer_list< float > lst)
    : InaVecALTIVEC(lst.begin()) {
    }

    inline void setFromScalar(const float val) {
        vec = vec_splats(val);
    }

    // Constructor from vec
    inline explicit InaVecALTIVEC(const float ptr[]) {
        vec = vec_xl(0, ptr);
    }

    inline InaVecALTIVEC& setFromArray(const float ptr[]) {
        vec = vec_xl(0, ptr);
        return *this;
    }

    inline InaVecALTIVEC& setFromAlignedArray(const float ptr[]) {
        vec = vec_ld(0, ptr);
        return *this;
    }

    inline InaVecALTIVEC& setFromIndirectArray(const float values[], const int inIndirection[]) {
        alignas(16) const std::array< float, 4 > tmp = { { values[inIndirection[0]],
                                                           values[inIndirection[1]],
                                                           values[inIndirection[2]],
                                                           values[inIndirection[3]] } };
        vec = vec_ld(0, &tmp[0]);
        return *this;
    }

    inline InaVecALTIVEC& setFromIndirect2DArray(const float inArray[], const int inIndirection1[],
                                                 const int inLeadingDimension, const int inIndirection2[]) {
        alignas(16) const std::array< float, 4 > tmp = { { inArray[inIndirection1[0] * inLeadingDimension + inIndirection2[0]],
                                                           inArray[inIndirection1[1] * inLeadingDimension + inIndirection2[1]],
                                                           inArray[inIndirection1[2] * inLeadingDimension + inIndirection2[2]],
                                                           inArray[inIndirection1[3] * inLeadingDimension + inIndirection2[3]] } };
        vec = vec_ld(0, &tmp[0]);
        return *this;
    }

    // Move back to array
    inline void storeInArray(float ptr[]) const {
        vec_xst(vec, 0, ptr);
    }

    inline void storeInAlignedArray(float ptr[]) const {
        vec_st(vec, 0, ptr);
    }

    // Acce to individual values
    inline float at(const int index) const {
        return vec_extract(vec, index);
    }

    // Horizontal operation
    inline float horizontalSum() const {
        __vector unsigned char perm2301 = { 0x8U, 0x9U, 0xAU, 0xBU, 0xCU, 0xDU, 0xEU, 0xFU,
                                            0x0U, 0x1U, 0x2U, 0x3U, 0x4U, 0x5U, 0x6U, 0x7U };
        __vector float middleres        = vec_add(vec, vec_perm(vec, vec, perm2301));
        __vector unsigned char perm1032 = { 0x04U, 0x05U, 0x06U, 0x07U, 0x00U, 0x01U, 0x02U, 0x03U,
                                            0x08U, 0x09U, 0x0AU, 0x0BU, 0x0CU, 0x0DU, 0x0EU, 0x0FU };
        __vector float res = vec_add(middleres, vec_perm(middleres, middleres, perm1032));
        return vec_extract(res, 0);
    }

    inline float horizontalMul() const {
        // Does vec_xxmadd could be faster
        __vector unsigned char perm2301 = { 0x8U, 0x9U, 0xAU, 0xBU, 0xCU, 0xDU, 0xEU, 0xFU,
                                            0x0U, 0x1U, 0x2U, 0x3U, 0x4U, 0x5U, 0x6U, 0x7U };
        __vector float middleres        = vec_mul(vec, vec_perm(vec, vec, perm2301));
        __vector unsigned char perm1032 = { 0x04U, 0x05U, 0x06U, 0x07U, 0x00U, 0x01U, 0x02U, 0x03U,
                                            0x08U, 0x09U, 0x0AU, 0x0BU, 0x0CU, 0x0DU, 0x0EU, 0x0FU };
        __vector float res = vec_mul(middleres, vec_perm(middleres, middleres, perm1032));
        return vec_extract(res, 0);
    }

    inline InaVecALTIVEC sqrt() const {
        return vec_sqrt(vec);
    }

    inline InaVecALTIVEC exp() const {
        const __vector float COEFF_LOG2E = vec_splats(float(InaFastExp::CoeffLog2E()));
        const __vector float COEFF_A     = vec_splats(float(InaFastExp::CoeffA32()));
        const __vector float COEFF_B     = vec_splats(float(InaFastExp::CoeffB32()));
        const __vector float COEFF_P5_A  = vec_splats(float(InaFastExp::GetCoefficient6_5()));
        const __vector float COEFF_P5_B  = vec_splats(float(InaFastExp::GetCoefficient6_4()));
        const __vector float COEFF_P5_C  = vec_splats(float(InaFastExp::GetCoefficient6_3()));
        const __vector float COEFF_P5_D  = vec_splats(float(InaFastExp::GetCoefficient6_2()));
        const __vector float COEFF_P5_E  = vec_splats(float(InaFastExp::GetCoefficient6_1()));
        const __vector float COEFF_P5_F  = vec_splats(float(InaFastExp::GetCoefficient6_0()));

        __vector float x = vec * COEFF_LOG2E;

        const __vector float fractional_part = x - InaVecALTIVEC(x).floor().vec;

        __vector float factor = (((((COEFF_P5_A * fractional_part + COEFF_P5_B) * fractional_part + COEFF_P5_C) * fractional_part + COEFF_P5_D) * fractional_part + COEFF_P5_E) * fractional_part + COEFF_P5_F);

        x -= factor;

        __vector int castedInteger = vec_cts(COEFF_A * x + COEFF_B, 0);

        return reinterpret_cast< __vector float >(castedInteger);
    }

    inline InaVecALTIVEC expLowAcc() const {
        const __vector float COEFF_LOG2E = vec_splats(float(InaFastExp::CoeffLog2E()));
        const __vector float COEFF_A     = vec_splats(float(InaFastExp::CoeffA32()));
        const __vector float COEFF_B     = vec_splats(float(InaFastExp::CoeffB32()));
        const __vector float COEFF_P5_D  = vec_splats(float(InaFastExp::GetCoefficient3_2()));
        const __vector float COEFF_P5_E  = vec_splats(float(InaFastExp::GetCoefficient3_1()));
        const __vector float COEFF_P5_F  = vec_splats(float(InaFastExp::GetCoefficient3_0()));

        __vector float x = vec * COEFF_LOG2E;

        const __vector float fractional_part = x - InaVecALTIVEC(x).floor().vec;

        __vector float factor = ((COEFF_P5_D * fractional_part + COEFF_P5_E) * fractional_part + COEFF_P5_F);

        x -= factor;

        __vector int castedInteger = vec_cts(COEFF_A * x + COEFF_B, 0);

        return reinterpret_cast< __vector float >(castedInteger);
    }

    inline InaVecALTIVEC rsqrt() const {
        return vec_rsqrt(vec);
    }

    inline InaVecALTIVEC abs() const {
        return vec_abs(vec);
    }

    inline InaVecALTIVEC floor() const {
        return vec_floor(vec);
    }

    inline InaVecALTIVEC signOf() const {
        const __vector float minus0 = reinterpret_cast< __vector float >(vec_splats(0x80000000U));
        const __vector float signs  = vec_and(vec, minus0);
        const __vector float gt0    = reinterpret_cast< __vector float >(vec_cmpeq(vec_splats(0.f), vec));
        return vec_and(vec_nand(gt0, gt0), vec_or(signs, vec_splats(1.f)));
    }

    inline InaVecALTIVEC isPositive() const {
        const __vector float testResult = reinterpret_cast< __vector float >(vec_cmpge(vec, vec_splats(0.f)));
        const __vector float ones       = vec_splats(1.f);
        return vec_and(testResult, ones);
    }

    inline InaVecALTIVEC isNegative() const {
        const __vector float testResult = reinterpret_cast< __vector float >(vec_cmpge(vec_splats(0.f), vec));
        const __vector float ones       = vec_splats(1.f);
        return vec_and(testResult, ones);
    }

    inline InaVecALTIVEC isPositiveStrict() const {
        const __vector float testResult = reinterpret_cast< __vector float >(vec_cmpgt(vec, vec_splats(0.f)));
        const __vector float ones       = vec_splats(1.f);
        return vec_and(testResult, ones);
    }

    inline InaVecALTIVEC isNegativeStrict() const {
        const __vector float testResult = reinterpret_cast< __vector float >(vec_cmpgt(vec_splats(0.f), vec));
        const __vector float ones       = vec_splats(1.f);
        return vec_and(testResult, ones);
    }

    inline InaVecALTIVEC isZero() const {
        const __vector float testResult = reinterpret_cast< __vector float >(vec_cmpeq(vec_splats(0.f), vec));
        const __vector float ones       = vec_splats(1.f);
        return vec_and(testResult, ones);
    }

    inline InaVecALTIVEC isNotZero() const {
        const __vector float testResult = reinterpret_cast< __vector float >(vec_cmpeq(vec_splats(0.f), vec));
        const __vector float ones       = vec_splats(1.f);
        return vec_and(vec_nand(testResult, testResult), ones);
    }

    inline InaVecMaskALTIVEC< float > isPositiveMask() const {
        return vec_cmpge(vec, vec_splats(0.f));
    }

    inline InaVecMaskALTIVEC< float > isNegativeMask() const {
        return vec_cmpge(vec_splats(0.f), vec);
    }

    inline InaVecMaskALTIVEC< float > isPositiveStrictMask() const {
        return vec_cmpgt(vec, vec_splats(0.f));
    }

    inline InaVecMaskALTIVEC< float > isNegativeStrictMask() const {
        return vec_cmpgt(vec_splats(0.f), vec);
        ;
    }

    inline InaVecMaskALTIVEC< float > isZeroMask() const {
        return vec_cmpeq(vec_splats(0.f), vec);
        ;
    }

    inline InaVecMaskALTIVEC< float > isNotZeroMask() const {
        return reinterpret_cast< __vector __bool int >(
            vec_nand(reinterpret_cast< __vector unsigned int >(vec_cmpeq(vec_splats(0.f), vec)),
                     reinterpret_cast< __vector unsigned int >(vec_splats(0xFFFFFFFFU))));
    }

    // Static basic methods
    inline static InaVecALTIVEC GetZero() {
        return InaVecALTIVEC(vec_splats(0.f));
    }

    inline static InaVecALTIVEC GetOne() {
        return InaVecALTIVEC(vec_splats(1.f));
    }

    inline static InaVecALTIVEC Min(const InaVecALTIVEC& inVec1, const InaVecALTIVEC& inVec2) {
        return vec_min(inVec1.vec, inVec2.vec);
    }

    inline static InaVecALTIVEC Max(const InaVecALTIVEC& inVec1, const InaVecALTIVEC& inVec2) {
        return vec_max(inVec1.vec, inVec2.vec);
    }

    inline static InaVecALTIVEC IsLowerOrEqual(const InaVecALTIVEC& inVec1, const InaVecALTIVEC& inVec2) {
        const __vector float testResult = reinterpret_cast< __vector float >(vec_cmpge(inVec2.vec, inVec1.vec));
        const __vector float ones       = vec_splats(1.f);
        return vec_and(testResult, ones);
    }

    inline static InaVecALTIVEC IsLower(const InaVecALTIVEC& inVec1, const InaVecALTIVEC& inVec2) {
        const __vector float testResult = reinterpret_cast< __vector float >(vec_cmpgt(inVec2.vec, inVec1.vec));
        const __vector float ones       = vec_splats(1.f);
        return vec_and(testResult, ones);
    }

    inline static InaVecALTIVEC IsGreaterOrEqual(const InaVecALTIVEC& inVec1, const InaVecALTIVEC& inVec2) {
        const __vector float testResult = reinterpret_cast< __vector float >(vec_cmpge(inVec1.vec, inVec2.vec));
        const __vector float ones       = vec_splats(1.f);
        return vec_and(testResult, ones);
    }

    inline static InaVecALTIVEC IsGreater(const InaVecALTIVEC& inVec1, const InaVecALTIVEC& inVec2) {
        const __vector float testResult = reinterpret_cast< __vector float >(vec_cmpgt(inVec1.vec, inVec2.vec));
        const __vector float ones       = vec_splats(1.f);
        return vec_and(testResult, ones);
    }

    inline static InaVecALTIVEC IsEqual(const InaVecALTIVEC& inVec1, const InaVecALTIVEC& inVec2) {
        const __vector float testResult = reinterpret_cast< __vector float >(vec_cmpeq(inVec1.vec, inVec2.vec));
        const __vector float ones       = vec_splats(1.f);
        return vec_and(testResult, ones);
    }

    inline static InaVecALTIVEC IsNotEqual(const InaVecALTIVEC& inVec1, const InaVecALTIVEC& inVec2) {
        const __vector float testResult = reinterpret_cast< __vector float >(vec_xor(reinterpret_cast< __vector unsigned >(vec_cmpeq(inVec1.vec, inVec2.vec)), vec_splats(0xFFFFFFFFU)));
        const __vector float ones       = vec_splats(1.f);
        return vec_and(testResult, ones);
    }

    inline static InaVecMaskALTIVEC< float > IsLowerOrEqualMask(const InaVecALTIVEC& inVec1, const InaVecALTIVEC& inVec2) {
        return vec_cmpge(inVec2.vec, inVec1.vec);
    }

    inline static InaVecMaskALTIVEC< float > IsLowerMask(const InaVecALTIVEC& inVec1, const InaVecALTIVEC& inVec2) {
        return vec_cmpgt(inVec2.vec, inVec1.vec);
    }

    inline static InaVecMaskALTIVEC< float > IsGreaterOrEqualMask(const InaVecALTIVEC& inVec1, const InaVecALTIVEC& inVec2) {
        return vec_cmpge(inVec1.vec, inVec2.vec);
    }

    inline static InaVecMaskALTIVEC< float > IsGreaterMask(const InaVecALTIVEC& inVec1, const InaVecALTIVEC& inVec2) {
        return vec_cmpgt(inVec1.vec, inVec2.vec);
    }

    inline static InaVecMaskALTIVEC< float > IsEqualMask(const InaVecALTIVEC& inVec1, const InaVecALTIVEC& inVec2) {
        return vec_cmpeq(inVec1.vec, inVec2.vec);
    }

    inline static InaVecMaskALTIVEC< float > IsNotEqualMask(const InaVecALTIVEC& inVec1, const InaVecALTIVEC& inVec2) {
        return vec_xor(reinterpret_cast< __vector __bool int >(vec_cmpeq(inVec1.vec, inVec2.vec)), reinterpret_cast< __vector __bool int >(vec_splats(0xFFFFFFFFU)));
    }

    inline static InaVecALTIVEC BitsAnd(const InaVecALTIVEC& inVec1, const InaVecALTIVEC& inVec2) {
        return vec_and(inVec1.vec, inVec2.vec);
    }

    inline static InaVecALTIVEC BitsNotAnd(const InaVecALTIVEC& inVec1, const InaVecALTIVEC& inVec2) {
        return vec_and(vec_nand(inVec1.vec, inVec1.vec), inVec2.vec);
    }

    inline static InaVecALTIVEC BitsOr(const InaVecALTIVEC& inVec1, const InaVecALTIVEC& inVec2) {
        return vec_or(inVec1.vec, inVec2.vec);
    }

    inline static InaVecALTIVEC BitsXor(const InaVecALTIVEC& inVec1, const InaVecALTIVEC& inVec2) {
        return vec_xor(inVec1.vec, inVec2.vec);
    }

    inline static const char* GetName() {
        return "InaVecALTIVEC<float>";
    }

    inline static InaIfElse< InaVecALTIVEC< float > >::ThenClass If(const InaVecMaskALTIVEC< float >& inTest) {
        return InaIfElse< InaVecALTIVEC< float > >::IfClass().If(inTest);
    }

    inline static InaVecALTIVEC IfElse(const InaVecMaskALTIVEC< float >& inMask, const InaVecALTIVEC& inIfTrue, const InaVecALTIVEC& inIfFalse) {
        return vec_or(IfTrue(inMask, inIfTrue.vec).vec,
                      IfFalse(inMask, inIfFalse.vec).vec);
    }

    inline static InaVecALTIVEC IfTrue(const InaVecMaskALTIVEC< float >& inMask, const InaVecALTIVEC& inIfTrue) {
        return vec_and(reinterpret_cast< __vector float >(inMask.getMask()), inIfTrue.vec);
    }

    inline static InaVecALTIVEC IfFalse(const InaVecMaskALTIVEC< float >& inMask, const InaVecALTIVEC& inIfFalse) {
        return reinterpret_cast< __vector float >(vec_and(vec_nand(reinterpret_cast< __vector unsigned int >(inMask.getMask()),
                                                                   reinterpret_cast< __vector unsigned int >(inMask.getMask())),
                                                          reinterpret_cast< __vector unsigned int >(inIfFalse.vec)));
    }

    // Inner operators
    inline InaVecALTIVEC< float >& operator+=(const InaVecALTIVEC< float >& inVec) {
        vec = vec_add(vec, inVec.vec);
        return *this;
    }

    inline InaVecALTIVEC< float >& operator-=(const InaVecALTIVEC< float >& inVec) {
        vec = vec_sub(vec, inVec.vec);
        return *this;
    }

    inline InaVecALTIVEC< float >& operator/=(const InaVecALTIVEC< float >& inVec) {
        vec = vec_div(vec, inVec.vec);
        return *this;
    }

    inline InaVecALTIVEC< float >& operator*=(const InaVecALTIVEC< float >& inVec) {
        vec = vec_mul(vec, inVec.vec);
        return *this;
    }

    inline InaVecALTIVEC< float > operator-() const {
        const __vector float minus0 = reinterpret_cast< __vector float >(vec_splats(0x80000000));
        return vec_xor(vec, minus0);
    }

    inline InaVecALTIVEC< float > pow(size_t power) const {
        return InaUtils::FastPow< InaVecALTIVEC< float > >(vec, power);
    }
};


// Bits operators
inline InaVecALTIVEC< float > operator&(const InaVecALTIVEC< float >& inVec1, const InaVecALTIVEC< float >& inVec2) {
    return InaVecALTIVEC< float >::BitsAnd(inVec1, inVec2);
}

inline InaVecALTIVEC< float > operator|(const InaVecALTIVEC< float >& inVec1, const InaVecALTIVEC< float >& inVec2) {
    return InaVecALTIVEC< float >::BitsOr(inVec1, inVec2);
}

inline InaVecALTIVEC< float > operator^(const InaVecALTIVEC< float >& inVec1, const InaVecALTIVEC< float >& inVec2) {
    return InaVecALTIVEC< float >::BitsXor(inVec1, inVec2);
}

// Dual operators
inline InaVecALTIVEC< float > operator+(const InaVecALTIVEC< float >& inVec1, const InaVecALTIVEC< float >& inVec2) {
    return vec_add(inVec1.getVec(), inVec2.getVec());
}

inline InaVecALTIVEC< float > operator-(const InaVecALTIVEC< float >& inVec1, const InaVecALTIVEC< float >& inVec2) {
    return vec_sub(inVec1.getVec(), inVec2.getVec());
}

inline InaVecALTIVEC< float > operator/(const InaVecALTIVEC< float >& inVec1, const InaVecALTIVEC< float >& inVec2) {
    return vec_div(inVec1.getVec(), inVec2.getVec());
}

inline InaVecALTIVEC< float > operator*(const InaVecALTIVEC< float >& inVec1, const InaVecALTIVEC< float >& inVec2) {
    return vec_mul(inVec1.getVec(), inVec2.getVec());
}

// Tests and comparions
inline InaVecMaskALTIVEC< float > operator<(const InaVecALTIVEC< float >& inVec1, const InaVecALTIVEC< float >& inVec2) {
    return InaVecALTIVEC< float >::IsLowerMask(inVec1, inVec2);
}

inline InaVecMaskALTIVEC< float > operator<=(const InaVecALTIVEC< float >& inVec1, const InaVecALTIVEC< float >& inVec2) {
    return InaVecALTIVEC< float >::IsLowerOrEqualMask(inVec1, inVec2);
}

inline InaVecMaskALTIVEC< float > operator>(const InaVecALTIVEC< float >& inVec1, const InaVecALTIVEC< float >& inVec2) {
    return InaVecALTIVEC< float >::IsGreaterMask(inVec1, inVec2);
}

inline InaVecMaskALTIVEC< float > operator>=(const InaVecALTIVEC< float >& inVec1, const InaVecALTIVEC< float >& inVec2) {
    return InaVecALTIVEC< float >::IsGreaterOrEqualMask(inVec1, inVec2);
}

inline InaVecMaskALTIVEC< float > operator==(const InaVecALTIVEC< float >& inVec1, const InaVecALTIVEC< float >& inVec2) {
    return InaVecALTIVEC< float >::IsEqualMask(inVec1, inVec2);
}

inline InaVecMaskALTIVEC< float > operator!=(const InaVecALTIVEC< float >& inVec1, const InaVecALTIVEC< float >& inVec2) {
    return InaVecALTIVEC< float >::IsNotEqualMask(inVec1, inVec2);
}


#endif
