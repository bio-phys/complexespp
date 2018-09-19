///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef INAVECSCALARDOUBLE_HPP
#define INAVECSCALARDOUBLE_HPP


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
class InaVecMaskSCALAR< double > {
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
inline InaVecMaskSCALAR< double > operator&(const InaVecMaskSCALAR< double >& inMask1, const InaVecMaskSCALAR< double >& inMask2) {
    return InaVecMaskSCALAR< double >::And(inMask1, inMask2);
}

inline InaVecMaskSCALAR< double > operator|(const InaVecMaskSCALAR< double >& inMask1, const InaVecMaskSCALAR< double >& inMask2) {
    return InaVecMaskSCALAR< double >::Or(inMask1, inMask2);
}

inline InaVecMaskSCALAR< double > operator^(const InaVecMaskSCALAR< double >& inMask1, const InaVecMaskSCALAR< double >& inMask2) {
    return InaVecMaskSCALAR< double >::Xor(inMask1, inMask2);
}

inline bool operator==(const InaVecMaskSCALAR< double >& inMask1, const InaVecMaskSCALAR< double >& inMask2) {
    return InaVecMaskSCALAR< double >::IsEqual(inMask1, inMask2);
}

inline bool operator!=(const InaVecMaskSCALAR< double >& inMask1, const InaVecMaskSCALAR< double >& inMask2) {
    return InaVecMaskSCALAR< double >::IsNotEqual(inMask1, inMask2);
}

// Vec type
template <>
class InaVecSCALAR< double > {
protected:
    double vec;

public:
    using VecRawType            = double;
    using MaskType              = InaVecMaskSCALAR< double >;
    using RealType              = double;
    static const int VecLength  = 1;
    static const int Alignement = 1;

    inline InaVecSCALAR() {
    }
    inline InaVecSCALAR(const InaVecSCALAR&) = default;
    inline InaVecSCALAR& operator=(const InaVecSCALAR&) = default;

    // Constructor from raw type
    inline /*not explicit*/ InaVecSCALAR(const double inVec)
    : vec(inVec) {
    }

    inline InaVecSCALAR& operator=(const double inVec) {
        vec = inVec;
        return *this;
    }

    inline void setFromRawType(const double inVec) {
        vec = inVec;
    }

    inline explicit operator double() const {
        return vec;
    }

    inline double getVec() const {
        return vec;
    }

    // Constructor from scalar

    // Constructor from vec
    inline InaVecSCALAR(const std::initializer_list< double > lst)
    : InaVecSCALAR(lst.begin()) {
    }

    inline explicit InaVecSCALAR(const double ptr[])
    : vec(*ptr) {
    }

    inline InaVecSCALAR& setFromArray(const double ptr[]) {
        vec = *ptr;
        return *this;
    }

    inline InaVecSCALAR setFromAlignedArray(const double ptr[]) {
        vec = *ptr;
        return *this;
    }

    inline InaVecSCALAR setFromIndirectArray(const double values[], const int inIndirection[]) {
        vec = values[inIndirection[0]];
        return *this;
    }

    inline InaVecSCALAR setFromIndirect2DArray(const double inArray[], const int inIndirection1[],
                                               const int inLeadingDimension, const int inIndirection2[]) {
        vec = inArray[inIndirection1[0] * inLeadingDimension + inIndirection2[0]];
        return *this;
    }

    // Move back to array
    inline void storeInArray(double ptr[]) const {
        ptr[0] = vec;
    }

    inline void storeInAlignedArray(double ptr[]) const {
        ptr[0] = vec;
    }

    // Acce to individual values
    inline double at(const int /*index*/) const {
        return vec;
    }

    // Horizontal operation
    inline double horizontalSum() const {
        return vec;
    }

    inline double horizontalMul() const {
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
        return vec < 0 ? -1 : vec > 0 ? 1 : 0;
    }

    inline InaVecSCALAR isPositive() const {
        return vec >= 0 ? 1 : 0;
    }

    inline InaVecSCALAR isNegative() const {
        return vec <= 0 ? 1 : 0;
    }

    inline InaVecSCALAR isPositiveStrict() const {
        return vec > 0 ? 1 : 0;
    }

    inline InaVecSCALAR isNegativeStrict() const {
        return vec < 0 ? 1 : 0;
    }

    inline InaVecSCALAR isZero() const {
        return vec == 0 ? 1 : 0;
    }

    inline InaVecSCALAR isNotZero() const {
        return vec == 0 ? 0 : 1;
    }

    inline InaVecMaskSCALAR< double > isPositiveMask() const {
        return vec >= 0 ? true : false;
    }

    inline InaVecMaskSCALAR< double > isNegativeMask() const {
        return vec <= 0 ? true : false;
    }

    inline InaVecMaskSCALAR< double > isPositiveStrictMask() const {
        return vec > 0 ? true : false;
    }

    inline InaVecMaskSCALAR< double > isNegativeStrictMask() const {
        return vec < 0 ? true : false;
    }

    inline InaVecMaskSCALAR< double > isZeroMask() const {
        return vec == 0 ? true : false;
    }

    inline InaVecMaskSCALAR< double > isNotZeroMask() const {
        return vec == 0 ? false : true;
    }

    // Static basic methods
    inline static InaVecSCALAR GetZero() {
        return 0;
    }

    inline static InaVecSCALAR GetOne() {
        return 1;
    }

    inline static InaVecSCALAR Min(const InaVecSCALAR& inVec1, const InaVecSCALAR& inVec2) {
        return (inVec1.vec <= inVec2.vec ? inVec1 : inVec2);
    }

    inline static InaVecSCALAR Max(const InaVecSCALAR& inVec1, const InaVecSCALAR& inVec2) {
        return (inVec1.vec >= inVec2.vec ? inVec1 : inVec2);
    }

    inline static InaVecSCALAR IsLowerOrEqual(const InaVecSCALAR& inVec1, const InaVecSCALAR& inVec2) {
        return inVec1.vec <= inVec2.vec ? 1 : 0;
    }

    inline static InaVecSCALAR IsLower(const InaVecSCALAR& inVec1, const InaVecSCALAR& inVec2) {
        return inVec1.vec < inVec2.vec ? 1 : 0;
    }

    inline static InaVecSCALAR IsGreaterOrEqual(const InaVecSCALAR& inVec1, const InaVecSCALAR& inVec2) {
        return inVec1.vec >= inVec2.vec ? 1 : 0;
    }

    inline static InaVecSCALAR IsGreater(const InaVecSCALAR& inVec1, const InaVecSCALAR& inVec2) {
        return inVec1.vec > inVec2.vec ? 1 : 0;
    }

    inline static InaVecSCALAR IsEqual(const InaVecSCALAR& inVec1, const InaVecSCALAR& inVec2) {
        return inVec1.vec == inVec2.vec ? 1 : 0;
    }

    inline static InaVecSCALAR IsNotEqual(const InaVecSCALAR& inVec1, const InaVecSCALAR& inVec2) {
        return inVec1.vec != inVec2.vec ? 1 : 0;
    }

    inline static InaVecMaskSCALAR< double > IsLowerOrEqualMask(const InaVecSCALAR& inVec1, const InaVecSCALAR& inVec2) {
        return inVec1.vec <= inVec2.vec ? true : false;
    }

    inline static InaVecMaskSCALAR< double > IsLowerMask(const InaVecSCALAR& inVec1, const InaVecSCALAR& inVec2) {
        return inVec1.vec < inVec2.vec ? true : false;
    }

    inline static InaVecMaskSCALAR< double > IsGreaterOrEqualMask(const InaVecSCALAR& inVec1, const InaVecSCALAR& inVec2) {
        return inVec1.vec >= inVec2.vec ? true : false;
    }

    inline static InaVecMaskSCALAR< double > IsGreaterMask(const InaVecSCALAR& inVec1, const InaVecSCALAR& inVec2) {
        return inVec1.vec > inVec2.vec ? true : false;
    }

    inline static InaVecMaskSCALAR< double > IsEqualMask(const InaVecSCALAR& inVec1, const InaVecSCALAR& inVec2) {
        return inVec1.vec == inVec2.vec ? true : false;
    }

    inline static InaVecMaskSCALAR< double > IsNotEqualMask(const InaVecSCALAR& inVec1, const InaVecSCALAR& inVec2) {
        return inVec1.vec != inVec2.vec ? true : false;
    }

    inline static InaVecSCALAR BitsAnd(const InaVecSCALAR& inVec1, const InaVecSCALAR& inVec2) {
        return InaUtils::ConvertBits< long int, double >(InaUtils::ConvertBits< double, long int >(inVec1.vec) & InaUtils::ConvertBits< double, long int >(inVec2.vec));
    }

    inline static InaVecSCALAR BitsNotAnd(const InaVecSCALAR& inVec1, const InaVecSCALAR& inVec2) {
        return InaUtils::ConvertBits< long int, double >((~InaUtils::ConvertBits< double, long int >(inVec1.vec)) & InaUtils::ConvertBits< double, long int >(inVec2.vec));
    }

    inline static InaVecSCALAR BitsOr(const InaVecSCALAR& inVec1, const InaVecSCALAR& inVec2) {
        return InaUtils::ConvertBits< long int, double >(InaUtils::ConvertBits< double, long int >(inVec1.vec) | InaUtils::ConvertBits< double, long int >(inVec2.vec));
    }

    inline static InaVecSCALAR BitsXor(const InaVecSCALAR& inVec1, const InaVecSCALAR& inVec2) {
        return InaUtils::ConvertBits< long int, double >(InaUtils::ConvertBits< double, long int >(inVec1.vec) ^ InaUtils::ConvertBits< double, long int >(inVec2.vec));
    }

    inline static const char* GetName() {
        return "InaVecSCALAR<double>";
    }

    inline static InaIfElse< InaVecSCALAR< double > >::ThenClass If(const InaVecMaskSCALAR< double >& inTest) {
        return InaIfElse< InaVecSCALAR< double > >::IfClass().If(inTest);
    }

    inline static InaVecSCALAR IfElse(const InaVecMaskSCALAR< double >& inMask, const InaVecSCALAR& inIfTrue, const InaVecSCALAR& inIfFalse) {
        return inMask ? inIfTrue : inIfFalse;
    }

    inline static InaVecSCALAR IfTrue(const InaVecMaskSCALAR< double >& inMask, const InaVecSCALAR& inIfTrue) {
        return inMask ? inIfTrue : 0;
    }

    inline static InaVecSCALAR IfFalse(const InaVecMaskSCALAR< double >& inMask, const InaVecSCALAR& inIfFalse) {
        return inMask ? 0 : inIfFalse;
    }

    // Inner operators
    inline InaVecSCALAR< double >& operator+=(const InaVecSCALAR< double >& inVec) {
        vec += inVec.vec;
        return *this;
    }

    inline InaVecSCALAR< double >& operator-=(const InaVecSCALAR< double >& inVec) {
        vec -= inVec.vec;
        return *this;
    }

    inline InaVecSCALAR< double >& operator/=(const InaVecSCALAR< double >& inVec) {
        vec /= inVec.vec;
        return *this;
    }

    inline InaVecSCALAR< double >& operator*=(const InaVecSCALAR< double >& inVec) {
        vec *= inVec.vec;
        return *this;
    }

    inline InaVecSCALAR< double > operator-() const {
        return -vec;
    }

    inline InaVecSCALAR< double > pow(size_t power) const {
        return InaUtils::FastPow< InaVecSCALAR< double > >(vec, power);
    }
};

// Bits operators
inline InaVecSCALAR< double > operator&(const InaVecSCALAR< double >& inVec1, const InaVecSCALAR< double >& inVec2) {
    return InaVecSCALAR< double >::BitsAnd(inVec1, inVec2);
}

inline InaVecSCALAR< double > operator|(const InaVecSCALAR< double >& inVec1, const InaVecSCALAR< double >& inVec2) {
    return InaVecSCALAR< double >::BitsOr(inVec1, inVec2);
}

inline InaVecSCALAR< double > operator^(const InaVecSCALAR< double >& inVec1, const InaVecSCALAR< double >& inVec2) {
    return InaVecSCALAR< double >::BitsXor(inVec1, inVec2);
}

// Dual operators
inline InaVecSCALAR< double > operator+(const InaVecSCALAR< double >& inVec1, const InaVecSCALAR< double >& inVec2) {
    return inVec1.getVec() + inVec2.getVec();
}

inline InaVecSCALAR< double > operator-(const InaVecSCALAR< double >& inVec1, const InaVecSCALAR< double >& inVec2) {
    return inVec1.getVec() - inVec2.getVec();
}

inline InaVecSCALAR< double > operator/(const InaVecSCALAR< double >& inVec1, const InaVecSCALAR< double >& inVec2) {
    return inVec1.getVec() / inVec2.getVec();
}

inline InaVecSCALAR< double > operator*(const InaVecSCALAR< double >& inVec1, const InaVecSCALAR< double >& inVec2) {
    return inVec1.getVec() * inVec2.getVec();
}

// Tests and comparions
inline InaVecMaskSCALAR< double > operator<(const InaVecSCALAR< double >& inVec1, const InaVecSCALAR< double >& inVec2) {
    return InaVecSCALAR< double >::IsLowerMask(inVec1, inVec2);
}

inline InaVecMaskSCALAR< double > operator<=(const InaVecSCALAR< double >& inVec1, const InaVecSCALAR< double >& inVec2) {
    return InaVecSCALAR< double >::IsLowerOrEqualMask(inVec1, inVec2);
}

inline InaVecMaskSCALAR< double > operator>(const InaVecSCALAR< double >& inVec1, const InaVecSCALAR< double >& inVec2) {
    return InaVecSCALAR< double >::IsGreaterMask(inVec1, inVec2);
}

inline InaVecMaskSCALAR< double > operator>=(const InaVecSCALAR< double >& inVec1, const InaVecSCALAR< double >& inVec2) {
    return InaVecSCALAR< double >::IsGreaterOrEqualMask(inVec1, inVec2);
}

inline InaVecMaskSCALAR< double > operator==(const InaVecSCALAR< double >& inVec1, const InaVecSCALAR< double >& inVec2) {
    return InaVecSCALAR< double >::IsEqualMask(inVec1, inVec2);
}

inline InaVecMaskSCALAR< double > operator!=(const InaVecSCALAR< double >& inVec1, const InaVecSCALAR< double >& inVec2) {
    return InaVecSCALAR< double >::IsNotEqualMask(inVec1, inVec2);
}


#endif // INAVECDOUBLE_HPP
