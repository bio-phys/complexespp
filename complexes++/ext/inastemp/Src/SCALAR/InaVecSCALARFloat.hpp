///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef INAVECSCALARFLOAT_HPP
#define INAVECSCALARFLOAT_HPP

#include "InastempConfig.h"
#include "Common/InaUtils.hpp"
#include "Common/InaIfElse.hpp"

#include <cmath>
#include <initializer_list>

template < class RealType >
class InaVecMaskSCALAR;

template < class RealType >
class InaVecSCALAR;

// Mask type
template <>
class InaVecMaskSCALAR< float > {
    bool mask;

public:
    // Classic constructors
    inline InaVecMaskSCALAR() {
    }

    inline InaVecMaskSCALAR(const InaVecMaskSCALAR&) = default;
    inline InaVecMaskSCALAR& operator=(const InaVecMaskSCALAR&) = default;

    // Native data type compatibility
    inline /*not explicit*/ InaVecMaskSCALAR(const bool inMask)
    : mask(inMask) {
    }

    inline InaVecMaskSCALAR& operator=(const bool inMask) {
        mask = inMask;
        return (*this);
    }

    inline explicit operator bool() const {
        return mask;
    }

    inline bool getMask() const {
        return mask;
    }

    // Binary methods
    inline InaVecMaskSCALAR Not() const {
        return InaVecMaskSCALAR(!mask);
    }

    inline bool isAllTrue() const {
        // true if all FF
        return mask;
    }

    inline bool isAllFalse() const {
        // true if all zero
        return !mask;
    }

    // Double args methods
    inline static InaVecMaskSCALAR And(const InaVecMaskSCALAR& inMask1, const InaVecMaskSCALAR& inMask2) {
        return InaVecMaskSCALAR(inMask1.mask & inMask2.mask);
    }

    inline static InaVecMaskSCALAR NotAnd(const InaVecMaskSCALAR& inMask1, const InaVecMaskSCALAR& inMask2) {
        return InaVecMaskSCALAR((!inMask1.mask) & inMask2.mask);
    }

    inline static InaVecMaskSCALAR Or(const InaVecMaskSCALAR& inMask1, const InaVecMaskSCALAR& inMask2) {
        return InaVecMaskSCALAR(inMask1.mask | inMask2.mask);
    }

    inline static InaVecMaskSCALAR Xor(const InaVecMaskSCALAR& inMask1, const InaVecMaskSCALAR& inMask2) {
        return InaVecMaskSCALAR(inMask1.mask ^ inMask2.mask);
    }

    inline static bool IsEqual(const InaVecMaskSCALAR& inMask1, const InaVecMaskSCALAR& inMask2) {
        return inMask1.mask == inMask2.mask;
    }

    inline static bool IsNotEqual(const InaVecMaskSCALAR& inMask1, const InaVecMaskSCALAR& inMask2) {
        return inMask1.mask != inMask2.mask;
    }
};

// Mask must have operators
inline InaVecMaskSCALAR< float > operator&(const InaVecMaskSCALAR< float >& inMask1, const InaVecMaskSCALAR< float >& inMask2) {
    return InaVecMaskSCALAR< float >::And(inMask1, inMask2);
}

inline InaVecMaskSCALAR< float > operator|(const InaVecMaskSCALAR< float >& inMask1, const InaVecMaskSCALAR< float >& inMask2) {
    return InaVecMaskSCALAR< float >::Or(inMask1, inMask2);
}

inline InaVecMaskSCALAR< float > operator^(const InaVecMaskSCALAR< float >& inMask1, const InaVecMaskSCALAR< float >& inMask2) {
    return InaVecMaskSCALAR< float >::Xor(inMask1, inMask2);
}

inline bool operator==(const InaVecMaskSCALAR< float >& inMask1, const InaVecMaskSCALAR< float >& inMask2) {
    return InaVecMaskSCALAR< float >::IsEqual(inMask1, inMask2);
}

inline bool operator!=(const InaVecMaskSCALAR< float >& inMask1, const InaVecMaskSCALAR< float >& inMask2) {
    return InaVecMaskSCALAR< float >::IsNotEqual(inMask1, inMask2);
}

// Vec type
template <>
class InaVecSCALAR< float > {
protected:
    float vec;

public:
    using VecRawType            = float;
    using MaskType              = InaVecMaskSCALAR< float >;
    using RealType              = float;
    static const int VecLength  = 1;
    static const int Alignement = 1;

    inline InaVecSCALAR() {
    }
    inline InaVecSCALAR(const InaVecSCALAR&) = default;
    inline InaVecSCALAR& operator=(const InaVecSCALAR&) = default;

    // Constructor from raw type
    inline /*not explicit*/ InaVecSCALAR(const float inVec)
    : vec(inVec) {
    }

    inline InaVecSCALAR& operator=(const float inVec) {
        vec = inVec;
        return *this;
    }

    inline InaVecSCALAR& setFromRawType(const float inVec) {
        vec = inVec;
        return *this;
    }

    inline explicit operator float() const {
        return vec;
    }

    inline float getVec() const {
        return vec;
    }

    // Constructor from scalar

    // Constructor from vec
    inline InaVecSCALAR(const std::initializer_list< float > lst)
    : InaVecSCALAR(lst.begin()) {
    }

    inline explicit InaVecSCALAR(const float ptr[])
    : vec(*ptr) {
    }

    inline InaVecSCALAR& setFromArray(const float ptr[]) {
        vec = *ptr;
        return *this;
    }

    inline InaVecSCALAR& setFromAlignedArray(const float ptr[]) {
        vec = *ptr;
        return *this;
    }

    inline InaVecSCALAR& setFromIndirectArray(const float values[], const int inIndirection[]) {
        vec = values[inIndirection[0]];
        return *this;
    }

    inline InaVecSCALAR& setFromIndirect2DArray(const float inArray[], const int inIndirection1[],
                                                const int inLeadingDimension, const int inIndirection2[]) {
        vec = inArray[inIndirection1[0] * inLeadingDimension + inIndirection2[0]];
        return *this;
    }

    // Move back to array
    inline void storeInArray(float ptr[]) const {
        ptr[0] = vec;
    }

    inline void storeInAlignedArray(float ptr[]) const {
        ptr[0] = vec;
    }

    // Acce to individual values
    inline float at(const int /*index*/) const {
        return vec;
    }

    // Horizontal operation
    inline float horizontalSum() const {
        return vec;
    }

    inline float horizontalMul() const {
        return vec;
    }

    inline InaVecSCALAR sqrt() const {
        return std::sqrt(vec);
    }

    inline InaVecSCALAR exp() const {
        return std::exp(vec);
    }

    inline InaVecSCALAR expLowAcc() const {
        return std::exp(vec);
    }

    inline InaVecSCALAR rsqrt() const {
        return 1 / std::sqrt(vec);
    }

    inline InaVecSCALAR abs() const {
        return vec < 0 ? -vec : vec;
    }

    inline InaVecSCALAR floor() const {
        return std::floor(vec);
    }

    inline InaVecSCALAR signOf() const {
        return vec < 0.f ? -1 : vec > 0.f ? 1.f : 0.f;
    }

    inline InaVecSCALAR isPositive() const {
        return vec >= 0.f ? 1.f : 0.f;
    }

    inline InaVecSCALAR isNegative() const {
        return vec <= 0.f ? 1.f : 0.f;
    }

    inline InaVecSCALAR isPositiveStrict() const {
        return vec > 0.f ? 1.f : 0.f;
    }

    inline InaVecSCALAR isNegativeStrict() const {
        return vec < 0.f ? 1.f : 0.f;
    }

    inline InaVecSCALAR isZero() const {
        return vec == 0.f ? 1.f : 0.f;
    }

    inline InaVecSCALAR isNotZero() const {
        return vec == 0.f ? 0.f : 1.f;
    }

    inline InaVecMaskSCALAR< float > isPositiveMask() const {
        return vec >= 0.f ? true : false;
    }

    inline InaVecMaskSCALAR< float > isNegativeMask() const {
        return vec <= 0.f ? true : false;
    }

    inline InaVecMaskSCALAR< float > isPositiveStrictMask() const {
        return vec > 0.f ? true : false;
    }

    inline InaVecMaskSCALAR< float > isNegativeStrictMask() const {
        return vec < 0.f ? true : false;
    }

    inline InaVecMaskSCALAR< float > isZeroMask() const {
        return vec == 0.f ? true : false;
    }

    inline InaVecMaskSCALAR< float > isNotZeroMask() const {
        return vec == 0.f ? false : true;
    }

    // Static basic methods
    inline static InaVecSCALAR GetZero() {
        return 0.f;
    }

    inline static InaVecSCALAR GetOne() {
        return 1.f;
    }

    inline static InaVecSCALAR Min(const InaVecSCALAR& inVec1, const InaVecSCALAR& inVec2) {
        return (inVec1.vec <= inVec2.vec ? inVec1 : inVec2);
    }

    inline static InaVecSCALAR Max(const InaVecSCALAR& inVec1, const InaVecSCALAR& inVec2) {
        return (inVec1.vec >= inVec2.vec ? inVec1 : inVec2);
    }

    inline static InaVecSCALAR IsLowerOrEqual(const InaVecSCALAR& inVec1, const InaVecSCALAR& inVec2) {
        return inVec1.vec <= inVec2.vec ? 1.f : 0.f;
    }

    inline static InaVecSCALAR IsLower(const InaVecSCALAR& inVec1, const InaVecSCALAR& inVec2) {
        return inVec1.vec < inVec2.vec ? 1.f : 0.f;
    }

    inline static InaVecSCALAR IsGreaterOrEqual(const InaVecSCALAR& inVec1, const InaVecSCALAR& inVec2) {
        return inVec1.vec >= inVec2.vec ? 1.f : 0.f;
    }

    inline static InaVecSCALAR IsGreater(const InaVecSCALAR& inVec1, const InaVecSCALAR& inVec2) {
        return inVec1.vec > inVec2.vec ? 1.f : 0.f;
    }

    inline static InaVecSCALAR IsEqual(const InaVecSCALAR& inVec1, const InaVecSCALAR& inVec2) {
        return inVec1.vec == inVec2.vec ? 1.f : 0.f;
    }

    inline static InaVecSCALAR IsNotEqual(const InaVecSCALAR& inVec1, const InaVecSCALAR& inVec2) {
        return inVec1.vec != inVec2.vec ? 1.f : 0.f;
    }

    inline static InaVecMaskSCALAR< float > IsLowerOrEqualMask(const InaVecSCALAR& inVec1, const InaVecSCALAR& inVec2) {
        return inVec1.vec <= inVec2.vec ? true : false;
    }

    inline static InaVecMaskSCALAR< float > IsLowerMask(const InaVecSCALAR& inVec1, const InaVecSCALAR& inVec2) {
        return inVec1.vec < inVec2.vec ? true : false;
    }

    inline static InaVecMaskSCALAR< float > IsGreaterOrEqualMask(const InaVecSCALAR& inVec1, const InaVecSCALAR& inVec2) {
        return inVec1.vec >= inVec2.vec ? true : false;
    }

    inline static InaVecMaskSCALAR< float > IsGreaterMask(const InaVecSCALAR& inVec1, const InaVecSCALAR& inVec2) {
        return inVec1.vec > inVec2.vec ? true : false;
    }

    inline static InaVecMaskSCALAR< float > IsEqualMask(const InaVecSCALAR& inVec1, const InaVecSCALAR& inVec2) {
        return inVec1.vec == inVec2.vec ? true : false;
    }

    inline static InaVecMaskSCALAR< float > IsNotEqualMask(const InaVecSCALAR& inVec1, const InaVecSCALAR& inVec2) {
        return inVec1.vec != inVec2.vec ? true : false;
    }

    inline static InaVecSCALAR BitsAnd(const InaVecSCALAR& inVec1, const InaVecSCALAR& inVec2) {
        return InaUtils::ConvertBits< int, float >(InaUtils::ConvertBits< float, int >(inVec1.vec) & InaUtils::ConvertBits< float, int >(inVec2.vec));
    }

    inline static InaVecSCALAR BitsNotAnd(const InaVecSCALAR& inVec1, const InaVecSCALAR& inVec2) {
        return InaUtils::ConvertBits< int, float >((~InaUtils::ConvertBits< float, int >(inVec1.vec)) & InaUtils::ConvertBits< float, int >(inVec2.vec));
    }

    inline static InaVecSCALAR BitsOr(const InaVecSCALAR& inVec1, const InaVecSCALAR& inVec2) {
        return InaUtils::ConvertBits< int, float >(InaUtils::ConvertBits< float, int >(inVec1.vec) | InaUtils::ConvertBits< float, int >(inVec2.vec));
    }

    inline static InaVecSCALAR BitsXor(const InaVecSCALAR& inVec1, const InaVecSCALAR& inVec2) {
        return InaUtils::ConvertBits< int, float >(InaUtils::ConvertBits< float, int >(inVec1.vec) ^ InaUtils::ConvertBits< float, int >(inVec2.vec));
    }

    inline static const char* GetName() {
        return "InaVecSCALAR<float>";
    }

    inline static InaIfElse< InaVecSCALAR< float > >::ThenClass If(const InaVecMaskSCALAR< float >& inTest) {
        return InaIfElse< InaVecSCALAR< float > >::IfClass().If(inTest);
    }

    inline static InaVecSCALAR IfElse(const InaVecMaskSCALAR< float >& inMask, const InaVecSCALAR& inIfTrue, const InaVecSCALAR& inIfFalse) {
        return inMask ? inIfTrue : inIfFalse;
    }

    inline static InaVecSCALAR IfTrue(const InaVecMaskSCALAR< float >& inMask, const InaVecSCALAR& inIfTrue) {
        return inMask ? inIfTrue : 0.f;
    }

    inline static InaVecSCALAR IfFalse(const InaVecMaskSCALAR< float >& inMask, const InaVecSCALAR& inIfFalse) {
        return inMask ? 0 : inIfFalse;
    }

    // Inner operators
    inline InaVecSCALAR< float >& operator+=(const InaVecSCALAR< float >& inVec) {
        vec += inVec.vec;
        return *this;
    }

    inline InaVecSCALAR< float >& operator-=(const InaVecSCALAR< float >& inVec) {
        vec -= inVec.vec;
        return *this;
    }

    inline InaVecSCALAR< float >& operator/=(const InaVecSCALAR< float >& inVec) {
        vec /= inVec.vec;
        return *this;
    }

    inline InaVecSCALAR< float >& operator*=(const InaVecSCALAR< float >& inVec) {
        vec *= inVec.vec;
        return *this;
    }

    inline InaVecSCALAR< float > operator-() const {
        return -vec;
    }

    inline InaVecSCALAR< float > pow(size_t power) const {
        return InaUtils::FastPow< InaVecSCALAR< float > >(vec, power);
    }
};

// Bits operators
inline InaVecSCALAR< float > operator&(const InaVecSCALAR< float >& inVec1, const InaVecSCALAR< float >& inVec2) {
    return InaVecSCALAR< float >::BitsAnd(inVec1, inVec2);
}

inline InaVecSCALAR< float > operator|(const InaVecSCALAR< float >& inVec1, const InaVecSCALAR< float >& inVec2) {
    return InaVecSCALAR< float >::BitsOr(inVec1, inVec2);
}

inline InaVecSCALAR< float > operator^(const InaVecSCALAR< float >& inVec1, const InaVecSCALAR< float >& inVec2) {
    return InaVecSCALAR< float >::BitsXor(inVec1, inVec2);
}

// Dual operators
inline InaVecSCALAR< float > operator+(const InaVecSCALAR< float >& inVec1, const InaVecSCALAR< float >& inVec2) {
    return inVec1.getVec() + inVec2.getVec();
}

inline InaVecSCALAR< float > operator-(const InaVecSCALAR< float >& inVec1, const InaVecSCALAR< float >& inVec2) {
    return inVec1.getVec() - inVec2.getVec();
}

inline InaVecSCALAR< float > operator/(const InaVecSCALAR< float >& inVec1, const InaVecSCALAR< float >& inVec2) {
    return inVec1.getVec() / inVec2.getVec();
}

inline InaVecSCALAR< float > operator*(const InaVecSCALAR< float >& inVec1, const InaVecSCALAR< float >& inVec2) {
    return inVec1.getVec() * inVec2.getVec();
}

// Tests and comparions
inline InaVecMaskSCALAR< float > operator<(const InaVecSCALAR< float >& inVec1, const InaVecSCALAR< float >& inVec2) {
    return InaVecSCALAR< float >::IsLowerMask(inVec1, inVec2);
}

inline InaVecMaskSCALAR< float > operator<=(const InaVecSCALAR< float >& inVec1, const InaVecSCALAR< float >& inVec2) {
    return InaVecSCALAR< float >::IsLowerOrEqualMask(inVec1, inVec2);
}

inline InaVecMaskSCALAR< float > operator>(const InaVecSCALAR< float >& inVec1, const InaVecSCALAR< float >& inVec2) {
    return InaVecSCALAR< float >::IsGreaterMask(inVec1, inVec2);
}

inline InaVecMaskSCALAR< float > operator>=(const InaVecSCALAR< float >& inVec1, const InaVecSCALAR< float >& inVec2) {
    return InaVecSCALAR< float >::IsGreaterOrEqualMask(inVec1, inVec2);
}

inline InaVecMaskSCALAR< float > operator==(const InaVecSCALAR< float >& inVec1, const InaVecSCALAR< float >& inVec2) {
    return InaVecSCALAR< float >::IsEqualMask(inVec1, inVec2);
}

inline InaVecMaskSCALAR< float > operator!=(const InaVecSCALAR< float >& inVec1, const InaVecSCALAR< float >& inVec2) {
    return InaVecSCALAR< float >::IsNotEqualMask(inVec1, inVec2);
}


#endif // INAVECFLOAT_HPP
