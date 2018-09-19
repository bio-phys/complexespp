///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef INAVECALTIVECDOUBLE_HPP
#define INAVECALTIVECDOUBLE_HPP

#include "InastempConfig.h"
#include "Common/InaIfElse.hpp"
#include "Common/InaUtils.hpp"

#ifndef INASTEMP_USE_ALTIVEC
#error InaVecALTIVEC<double> is included but ALTIVEC is not enable in the configuration
#endif

#include "Common/InaFastExp.hpp"

#include <altivec.h>
#undef bool
#undef vector
#undef pixel

#include <cmath>
#include <initializer_list>

// Forward declarations
template < class RealType >
class InaVecMaskALTIVEC;

template < class RealType >
class InaVecALTIVEC;

// Mask type
template <>
class alignas(16) InaVecMaskALTIVEC< double > {
    __vector __bool long long mask;

public:
    // Classic constructors
    inline InaVecMaskALTIVEC() {
    }

    InaVecMaskALTIVEC(const InaVecMaskALTIVEC&) = default;
    inline InaVecMaskALTIVEC& operator=(const InaVecMaskALTIVEC&) = default;

    // Native data type compatibility
    inline /*not explicit*/ InaVecMaskALTIVEC(const __vector __bool long long inMask)
    : mask(inMask) {
    }

    inline InaVecMaskALTIVEC& operator=(const __vector __bool long long inMask) {
        mask = inMask;
        return (*this);
    }

    inline explicit operator __vector __bool long long() const {
        return mask;
    }

    inline __vector __bool long long getMask() const {
        return mask;
    }

    // Bool data type compatibility
    inline explicit InaVecMaskALTIVEC(const bool inBool) {
        const __vector __bool long long tmpMaskFF = reinterpret_cast< __vector __bool long long >(vec_splats(0xFFFFFFFFFFFFFFFFULL));
        mask                                      = (inBool ? tmpMaskFF : reinterpret_cast< __vector __bool long long >(vec_splats(0x0ULL)));
    }

    inline InaVecMaskALTIVEC& operator=(const bool inBool) {
        const __vector __bool long long tmpMaskFF = reinterpret_cast< __vector __bool long long >(vec_splats(0xFFFFFFFFFFFFFFFFULL));
        mask                                      = (inBool ? tmpMaskFF : reinterpret_cast< __vector __bool long long >(vec_splats(0x0ULL)));
        return (*this);
    }

    // Binary methods
    inline InaVecMaskALTIVEC Not() const {
        const __vector __bool long long tmpMaskFF = reinterpret_cast< __vector __bool long long >(vec_splats(0xFFFFFFFFFFFFFFFFULL));
        return NotAnd(mask, tmpMaskFF);
    }

    inline bool isAllTrue() const {
        // true if all FF => !FF => 0 & FF => 0
        const __vector __bool long long tmpMaskFF = reinterpret_cast< __vector __bool long long >(vec_splats(0xFFFFFFFFFFFFFFFFULL));
        ;
        const int res = vec_all_eq(mask, tmpMaskFF);
        return static_cast< bool >(res);
    }

    inline bool isAllFalse() const {
        // true if all zero
        const int res = vec_all_eq(mask, reinterpret_cast< __vector __bool long long >(vec_xor(reinterpret_cast< __vector unsigned int >(mask),
                                                                                               reinterpret_cast< __vector unsigned int >(mask))));
        return static_cast< bool >(res);
    }

    // Double args methods
    inline static InaVecMaskALTIVEC And(const InaVecMaskALTIVEC& inMask1, const InaVecMaskALTIVEC& inMask2) {
        return InaVecMaskALTIVEC(reinterpret_cast< __vector __bool long long >(
            vec_and(reinterpret_cast< __vector unsigned char >(inMask1.mask),
                    reinterpret_cast< __vector unsigned char >(inMask2.mask))));
    }

    inline static InaVecMaskALTIVEC NotAnd(const InaVecMaskALTIVEC& inMask1, const InaVecMaskALTIVEC& inMask2) {
        return InaVecMaskALTIVEC(reinterpret_cast< __vector __bool long long >(
            vec_nand(reinterpret_cast< __vector unsigned char >(inMask1.mask),
                     reinterpret_cast< __vector unsigned char >(inMask2.mask))));
    }

    inline static InaVecMaskALTIVEC Or(const InaVecMaskALTIVEC& inMask1, const InaVecMaskALTIVEC& inMask2) {
        return InaVecMaskALTIVEC(reinterpret_cast< __vector __bool long long >(
            vec_or(reinterpret_cast< __vector unsigned char >(inMask1.mask),
                   reinterpret_cast< __vector unsigned char >(inMask2.mask))));
    }

    inline static InaVecMaskALTIVEC Xor(const InaVecMaskALTIVEC& inMask1, const InaVecMaskALTIVEC& inMask2) {
        return InaVecMaskALTIVEC(reinterpret_cast< __vector __bool long long >(
            vec_xor(reinterpret_cast< __vector unsigned char >(inMask1.mask),
                    reinterpret_cast< __vector unsigned char >(inMask2.mask))));
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
inline InaVecMaskALTIVEC< double > operator&(const InaVecMaskALTIVEC< double >& inMask1, const InaVecMaskALTIVEC< double >& inMask2) {
    return InaVecMaskALTIVEC< double >::And(inMask1, inMask2);
}

inline InaVecMaskALTIVEC< double > operator|(const InaVecMaskALTIVEC< double >& inMask1, const InaVecMaskALTIVEC< double >& inMask2) {
    return InaVecMaskALTIVEC< double >::Or(inMask1, inMask2);
}

inline InaVecMaskALTIVEC< double > operator^(const InaVecMaskALTIVEC< double >& inMask1, const InaVecMaskALTIVEC< double >& inMask2) {
    return InaVecMaskALTIVEC< double >::Xor(inMask1, inMask2);
}

inline bool operator==(const InaVecMaskALTIVEC< double >& inMask1, const InaVecMaskALTIVEC< double >& inMask2) {
    return InaVecMaskALTIVEC< double >::IsEqual(inMask1, inMask2);
}

inline bool operator!=(const InaVecMaskALTIVEC< double >& inMask1, const InaVecMaskALTIVEC< double >& inMask2) {
    return InaVecMaskALTIVEC< double >::IsNotEqual(inMask1, inMask2);
}

// Vec type
template <>
class alignas(16) InaVecALTIVEC< double > {
protected:
    __vector double vec;

public:
    using VecRawType            = __vector double;
    using MaskType              = InaVecMaskALTIVEC< double >;
    using RealType              = double;
    static const int VecLength  = 2;
    static const int Alignement = 16;

    inline InaVecALTIVEC() {
    }
    inline InaVecALTIVEC(const InaVecALTIVEC&) = default;
    inline InaVecALTIVEC& operator=(const InaVecALTIVEC&) = default;

    // Constructor from raw type
    inline /*not explicit*/ InaVecALTIVEC(const __vector double inVec)
    : vec(inVec) {
    }

    inline InaVecALTIVEC& operator=(const __vector double inVec) {
        vec = inVec;
        return *this;
    }

    inline void setFromRawType(const __vector double inVec) {
        vec = inVec;
    }

    inline explicit operator __vector double() const {
        return vec;
    }

    inline __vector double getVec() const {
        return vec;
    }

    // Constructor from scalar
    inline /*not explicit*/ InaVecALTIVEC(const double val)
    : vec(vec_splats(val)) {
    }

    inline InaVecALTIVEC& operator=(const double val) {
        vec = vec_splats(val);
        return *this;
    }

    // Constructor from vec
    inline InaVecALTIVEC(const std::initializer_list< double > lst)
    : InaVecALTIVEC(lst.begin()) {
    }

    inline void setFromScalar(const double val) {
        vec = vec_splats(val);
    }

    // Constructor from vec
    inline explicit InaVecALTIVEC(const double ptr[]) {
        vec = vec_xl(0, const_cast< double* >(ptr));
    }

    inline InaVecALTIVEC& setFromArray(const double ptr[]) {
        vec = vec_xl(0, const_cast< double* >(ptr));
        return *this;
    }

    inline InaVecALTIVEC& setFromAlignedArray(const double ptr[]) {
        // TODO use vec_ld without cast?
        vec = reinterpret_cast< __vector double >(vec_ld(0, reinterpret_cast< const unsigned long* >(&ptr[0])));
        return *this;
    }

    inline InaVecALTIVEC& setFromIndirectArray(const double values[], const int inIndirection[]) {
        alignas(16) const std::array< double, 2 > tmp = { { values[inIndirection[0]],
                                                            values[inIndirection[1]] } };
        vec = reinterpret_cast< __vector double >(vec_ld(0, reinterpret_cast< const unsigned long* >(&tmp[0])));
        return *this;
    }

    inline InaVecALTIVEC& setFromIndirect2DArray(const double inArray[], const int inIndirection1[],
                                                 const int inLeadingDimension, const int inIndirection2[]) {
        alignas(16) const std::array< double, 2 > tmp = { { inArray[inIndirection1[0] * inLeadingDimension + inIndirection2[0]],
                                                            inArray[inIndirection1[1] * inLeadingDimension + inIndirection2[1]] } };
        // TODO use vec_ld without cast?
        vec = reinterpret_cast< __vector double >(vec_ld(0, reinterpret_cast< const unsigned long* >(&tmp[0])));
        return *this;
    }

    // Move back to array
    inline void storeInArray(double ptr[]) const {
        // TODO remove cast
        vec_xst(reinterpret_cast< __vector unsigned int >(vec), 0, reinterpret_cast< unsigned int* >(ptr));
    }

    inline void storeInAlignedArray(double ptr[]) const {
        // TODO remove cast
        vec_st(reinterpret_cast< __vector unsigned int >(vec), 0, reinterpret_cast< unsigned int* >(ptr));
    }

    // Acce to individual values
    inline double at(const int index) const {
        return vec_extract(vec, index);
    }

    // Horizontal operation
    inline double horizontalSum() const {
        __vector unsigned char perm2301 = { 0x8U, 0x9U, 0xAU, 0xBU, 0xCU, 0xDU, 0xEU, 0xFU,
                                            0x0U, 0x1U, 0x2U, 0x3U, 0x4U, 0x5U, 0x6U, 0x7U };
        __vector double res = vec_add(vec, vec_perm(vec, vec, perm2301));
        return vec_extract(res, 0);
    }

    inline double horizontalMul() const {
        // Does vec_xxmadd could be faster
        __vector unsigned char perm2301 = { 0x8U, 0x9U, 0xAU, 0xBU, 0xCU, 0xDU, 0xEU, 0xFU,
                                            0x0U, 0x1U, 0x2U, 0x3U, 0x4U, 0x5U, 0x6U, 0x7U };
        __vector double res = vec_mul(vec, vec_perm(vec, vec, perm2301));
        return vec_extract(res, 0);
    }

    inline InaVecALTIVEC sqrt() const {
        return vec_sqrt(vec);
    }

    inline InaVecALTIVEC exp() const {
        const __vector double COEFF_LOG2E = vec_splats(double(InaFastExp::CoeffLog2E()));
        const __vector double COEFF_A     = vec_splats(double(InaFastExp::CoeffA64()));
        const __vector double COEFF_B     = vec_splats(double(InaFastExp::CoeffB64()));
        const __vector double COEFF_P5_X  = vec_splats(double(InaFastExp::GetCoefficient9_8()));
        const __vector double COEFF_P5_Y  = vec_splats(double(InaFastExp::GetCoefficient9_7()));
        const __vector double COEFF_P5_Z  = vec_splats(double(InaFastExp::GetCoefficient9_6()));
        const __vector double COEFF_P5_A  = vec_splats(double(InaFastExp::GetCoefficient9_5()));
        const __vector double COEFF_P5_B  = vec_splats(double(InaFastExp::GetCoefficient9_4()));
        const __vector double COEFF_P5_C  = vec_splats(double(InaFastExp::GetCoefficient9_3()));
        const __vector double COEFF_P5_D  = vec_splats(double(InaFastExp::GetCoefficient9_2()));
        const __vector double COEFF_P5_E  = vec_splats(double(InaFastExp::GetCoefficient9_1()));
        const __vector double COEFF_P5_F  = vec_splats(double(InaFastExp::GetCoefficient9_0()));

        __vector double x = vec * COEFF_LOG2E;

        const __vector double fractional_part = x - InaVecALTIVEC(x).floor().vec;

        __vector double factor = ((((((((COEFF_P5_X * fractional_part + COEFF_P5_Y) * fractional_part + COEFF_P5_Z) * fractional_part + COEFF_P5_A) * fractional_part + COEFF_P5_B) * fractional_part + COEFF_P5_C) * fractional_part + COEFF_P5_D) * fractional_part + COEFF_P5_E) * fractional_part + COEFF_P5_F);

        x -= factor;

        x = COEFF_A * x + COEFF_B;

        // TODO find conversion function
        //__vector long castedInteger = vec_ctsl(x, 0);
        //return reinterpret_cast<__vector double>(castedInteger);
        alignas(16) double tmpptr[2];
        vec_st(reinterpret_cast< __vector unsigned int >(x), 0, reinterpret_cast< unsigned int* >(tmpptr));

#ifdef INASTEMP_USE_XL
        alignas(16) long long ltmpptr[2];
#else
        alignas(16) long ltmpptr[2];
#endif
        ltmpptr[0] = long(tmpptr[0]);
        ltmpptr[1] = long(tmpptr[1]);
        return reinterpret_cast< __vector double >(vec_xl(0, ltmpptr));
    }

    inline InaVecALTIVEC expLowAcc() const {
        const __vector double COEFF_LOG2E = vec_splats(double(InaFastExp::CoeffLog2E()));
        const __vector double COEFF_A     = vec_splats(double(InaFastExp::CoeffA64()));
        const __vector double COEFF_B     = vec_splats(double(InaFastExp::CoeffB64()));
        const __vector double COEFF_P5_C  = vec_splats(double(InaFastExp::GetCoefficient4_3()));
        const __vector double COEFF_P5_D  = vec_splats(double(InaFastExp::GetCoefficient4_2()));
        const __vector double COEFF_P5_E  = vec_splats(double(InaFastExp::GetCoefficient4_1()));
        const __vector double COEFF_P5_F  = vec_splats(double(InaFastExp::GetCoefficient4_0()));

        __vector double x = vec * COEFF_LOG2E;

        const __vector double fractional_part = x - InaVecALTIVEC(x).floor().vec;

        __vector double factor = (((COEFF_P5_C * fractional_part + COEFF_P5_D) * fractional_part + COEFF_P5_E) * fractional_part + COEFF_P5_F);

        x -= factor;

        x = COEFF_A * x + COEFF_B;

        // TODO find conversion function
        //__vector long castedInteger = vec_cts(x, 0);
        //return reinterpret_cast<__vector double>(castedInteger);
        alignas(16) double tmpptr[2];
        vec_st(reinterpret_cast< __vector unsigned int >(x), 0, reinterpret_cast< unsigned int* >(tmpptr));

#ifdef INASTEMP_USE_XL
        alignas(16) long long ltmpptr[2];
#else
        alignas(16) long ltmpptr[2];
#endif
        ltmpptr[0] = long(tmpptr[0]);
        ltmpptr[1] = long(tmpptr[1]);
        return reinterpret_cast< __vector double >(vec_xl(0, ltmpptr));
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
        const __vector double minus0 = reinterpret_cast< __vector double >(vec_splats(0x8000000000000000ULL));
        const __vector double signs  = vec_and(vec, minus0);
        const __vector double ge0    = reinterpret_cast< __vector double >(vec_cmpeq(vec_splats(0.), vec));
        return vec_and(vec_nand(ge0, ge0), vec_or(signs, vec_splats(1.)));
    }

    inline InaVecALTIVEC isPositive() const {
        const __vector double testResult = reinterpret_cast< __vector double >(vec_cmpge(vec, vec_splats(0.)));
        const __vector double ones       = vec_splats(1.);
        return vec_and(testResult, ones);
    }

    inline InaVecALTIVEC isNegative() const {
        const __vector double testResult = reinterpret_cast< __vector double >(vec_cmpge(vec_splats(0.), vec));
        const __vector double ones       = vec_splats(1.);
        return vec_and(testResult, ones);
    }

    inline InaVecALTIVEC isPositiveStrict() const {
        const __vector double testResult = reinterpret_cast< __vector double >(vec_cmpgt(vec, vec_splats(0.)));
        const __vector double ones       = vec_splats(1.);
        return vec_and(testResult, ones);
    }

    inline InaVecALTIVEC isNegativeStrict() const {
        const __vector double testResult = reinterpret_cast< __vector double >(vec_cmpgt(vec_splats(0.), vec));
        const __vector double ones       = vec_splats(1.);
        return vec_and(testResult, ones);
    }

    inline InaVecALTIVEC isZero() const {
        const __vector double testResult = reinterpret_cast< __vector double >(vec_cmpeq(vec_splats(0.), vec));
        const __vector double ones       = vec_splats(1.);
        return vec_and(testResult, ones);
    }

    inline InaVecALTIVEC isNotZero() const {
        const __vector double testResult = reinterpret_cast< __vector double >(vec_cmpeq(vec_splats(0.), vec));
        const __vector double ones       = vec_splats(1.);
        return vec_and(vec_nand(testResult, testResult), ones);
    }

    inline InaVecMaskALTIVEC< double > isPositiveMask() const {
        return vec_cmpge(vec, vec_splats(0.));
    }

    inline InaVecMaskALTIVEC< double > isNegativeMask() const {
        return vec_cmpge(vec_splats(0.), vec);
    }

    inline InaVecMaskALTIVEC< double > isPositiveStrictMask() const {
        return vec_cmpgt(vec, vec_splats(0.));
    }

    inline InaVecMaskALTIVEC< double > isNegativeStrictMask() const {
        return vec_cmpgt(vec_splats(0.), vec);
        ;
    }

    inline InaVecMaskALTIVEC< double > isZeroMask() const {
        return vec_cmpeq(vec_splats(0.), vec);
        ;
    }

    inline InaVecMaskALTIVEC< double > isNotZeroMask() const {
        return reinterpret_cast< __vector __bool long long >(
            vec_nand(reinterpret_cast< __vector unsigned char >(vec_cmpeq(vec_splats(0.), vec)),
                     reinterpret_cast< __vector unsigned char >(vec_splats(0xFFFFFFFFU))));
    }

    // Static basic methods
    inline static InaVecALTIVEC GetZero() {
        return InaVecALTIVEC(vec_splats(0.));
    }

    inline static InaVecALTIVEC GetOne() {
        return InaVecALTIVEC(vec_splats(1.));
    }

    inline static InaVecALTIVEC Min(const InaVecALTIVEC& inVec1, const InaVecALTIVEC& inVec2) {
        return vec_min(inVec1.vec, inVec2.vec);
    }

    inline static InaVecALTIVEC Max(const InaVecALTIVEC& inVec1, const InaVecALTIVEC& inVec2) {
        return vec_max(inVec1.vec, inVec2.vec);
    }

    inline static InaVecALTIVEC IsLowerOrEqual(const InaVecALTIVEC& inVec1, const InaVecALTIVEC& inVec2) {
        const __vector double testResult = reinterpret_cast< __vector double >(vec_cmpge(inVec2.vec, inVec1.vec));
        const __vector double ones       = vec_splats(1.);
        return vec_and(testResult, ones);
    }

    inline static InaVecALTIVEC IsLower(const InaVecALTIVEC& inVec1, const InaVecALTIVEC& inVec2) {
        const __vector double testResult = reinterpret_cast< __vector double >(vec_cmpgt(inVec2.vec, inVec1.vec));
        const __vector double ones       = vec_splats(1.);
        return vec_and(testResult, ones);
    }

    inline static InaVecALTIVEC IsGreaterOrEqual(const InaVecALTIVEC& inVec1, const InaVecALTIVEC& inVec2) {
        const __vector double testResult = reinterpret_cast< __vector double >(vec_cmpge(inVec1.vec, inVec2.vec));
        const __vector double ones       = vec_splats(1.);
        return vec_and(testResult, ones);
    }

    inline static InaVecALTIVEC IsGreater(const InaVecALTIVEC& inVec1, const InaVecALTIVEC& inVec2) {
        const __vector double testResult = reinterpret_cast< __vector double >(vec_cmpgt(inVec1.vec, inVec2.vec));
        const __vector double ones       = vec_splats(1.);
        return vec_and(testResult, ones);
    }

    inline static InaVecALTIVEC IsEqual(const InaVecALTIVEC& inVec1, const InaVecALTIVEC& inVec2) {
        const __vector double testResult = reinterpret_cast< __vector double >(vec_cmpeq(inVec1.vec, inVec2.vec));
        const __vector double ones       = vec_splats(1.);
        return vec_and(testResult, ones);
    }

    inline static InaVecALTIVEC IsNotEqual(const InaVecALTIVEC& inVec1, const InaVecALTIVEC& inVec2) {
        const __vector double testResult = reinterpret_cast< __vector double >(vec_xor(reinterpret_cast< __vector unsigned >(vec_cmpeq(inVec1.vec, inVec2.vec)), vec_splats(0xFFFFFFFFU)));
        const __vector double ones       = vec_splats(1.);
        return vec_and(testResult, ones);
    }

    inline static InaVecMaskALTIVEC< double > IsLowerOrEqualMask(const InaVecALTIVEC& inVec1, const InaVecALTIVEC& inVec2) {
        return vec_cmpge(inVec2.vec, inVec1.vec);
    }

    inline static InaVecMaskALTIVEC< double > IsLowerMask(const InaVecALTIVEC& inVec1, const InaVecALTIVEC& inVec2) {
        return vec_cmpgt(inVec2.vec, inVec1.vec);
    }

    inline static InaVecMaskALTIVEC< double > IsGreaterOrEqualMask(const InaVecALTIVEC& inVec1, const InaVecALTIVEC& inVec2) {
        return vec_cmpge(inVec1.vec, inVec2.vec);
    }

    inline static InaVecMaskALTIVEC< double > IsGreaterMask(const InaVecALTIVEC& inVec1, const InaVecALTIVEC& inVec2) {
        return vec_cmpgt(inVec1.vec, inVec2.vec);
    }

    inline static InaVecMaskALTIVEC< double > IsEqualMask(const InaVecALTIVEC& inVec1, const InaVecALTIVEC& inVec2) {
        return vec_cmpeq(inVec1.vec, inVec2.vec);
    }

    inline static InaVecMaskALTIVEC< double > IsNotEqualMask(const InaVecALTIVEC& inVec1, const InaVecALTIVEC& inVec2) {
        return reinterpret_cast< __vector __bool long long >(vec_xor(reinterpret_cast< __vector unsigned int >(vec_cmpeq(inVec1.vec, inVec2.vec)),
                                                                     reinterpret_cast< __vector unsigned int >(vec_splats(0xFFFFFFFFU))));
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
        return "InaVecALTIVEC<double>";
    }

    inline static InaIfElse< InaVecALTIVEC< double > >::ThenClass If(const InaVecMaskALTIVEC< double >& inTest) {
        return InaIfElse< InaVecALTIVEC< double > >::IfClass().If(inTest);
    }

    inline static InaVecALTIVEC IfElse(const InaVecMaskALTIVEC< double >& inMask, const InaVecALTIVEC& inIfTrue, const InaVecALTIVEC& inIfFalse) {
        return vec_or(IfTrue(inMask, inIfTrue.vec).vec,
                      IfFalse(inMask, inIfFalse.vec).vec);
    }

    inline static InaVecALTIVEC IfTrue(const InaVecMaskALTIVEC< double >& inMask, const InaVecALTIVEC& inIfTrue) {
        return vec_and(reinterpret_cast< __vector double >(inMask.getMask()), inIfTrue.vec);
    }

    inline static InaVecALTIVEC IfFalse(const InaVecMaskALTIVEC< double >& inMask, const InaVecALTIVEC& inIfFalse) {
        return reinterpret_cast< __vector double >(vec_and(vec_nand(reinterpret_cast< __vector unsigned int >(inMask.getMask()),
                                                                    reinterpret_cast< __vector unsigned int >(inMask.getMask())),
                                                           reinterpret_cast< __vector unsigned int >(inIfFalse.vec)));
    }

    // Inner operators
    inline InaVecALTIVEC< double >& operator+=(const InaVecALTIVEC< double >& inVec) {
        vec = vec_add(vec, inVec.vec);
        return *this;
    }

    inline InaVecALTIVEC< double >& operator-=(const InaVecALTIVEC< double >& inVec) {
        vec = vec_sub(vec, inVec.vec);
        return *this;
    }

    inline InaVecALTIVEC< double >& operator/=(const InaVecALTIVEC< double >& inVec) {
        vec = vec_div(vec, inVec.vec);
        return *this;
    }

    inline InaVecALTIVEC< double >& operator*=(const InaVecALTIVEC< double >& inVec) {
        vec = vec_mul(vec, inVec.vec);
        return *this;
    }

    inline InaVecALTIVEC< double > operator-() const {
        const __vector double minus0 = reinterpret_cast< __vector double >(vec_splats(0x8000000000000000ULL));
        return vec_xor(vec, minus0);
    }

    inline InaVecALTIVEC< double > pow(size_t power) const {
        return InaUtils::FastPow< InaVecALTIVEC< double > >(vec, power);
    }
};

// Bits operators
inline InaVecALTIVEC< double > operator&(const InaVecALTIVEC< double >& inVec1, const InaVecALTIVEC< double >& inVec2) {
    return InaVecALTIVEC< double >::BitsAnd(inVec1, inVec2);
}

inline InaVecALTIVEC< double > operator|(const InaVecALTIVEC< double >& inVec1, const InaVecALTIVEC< double >& inVec2) {
    return InaVecALTIVEC< double >::BitsOr(inVec1, inVec2);
}

inline InaVecALTIVEC< double > operator^(const InaVecALTIVEC< double >& inVec1, const InaVecALTIVEC< double >& inVec2) {
    return InaVecALTIVEC< double >::BitsXor(inVec1, inVec2);
}

// Dual operators
inline InaVecALTIVEC< double > operator+(const InaVecALTIVEC< double >& inVec1, const InaVecALTIVEC< double >& inVec2) {
    return vec_add(inVec1.getVec(), inVec2.getVec());
}

inline InaVecALTIVEC< double > operator-(const InaVecALTIVEC< double >& inVec1, const InaVecALTIVEC< double >& inVec2) {
    return vec_sub(inVec1.getVec(), inVec2.getVec());
}

inline InaVecALTIVEC< double > operator/(const InaVecALTIVEC< double >& inVec1, const InaVecALTIVEC< double >& inVec2) {
    return vec_div(inVec1.getVec(), inVec2.getVec());
}

inline InaVecALTIVEC< double > operator*(const InaVecALTIVEC< double >& inVec1, const InaVecALTIVEC< double >& inVec2) {
    return vec_mul(inVec1.getVec(), inVec2.getVec());
}

// Tests and comparions
inline InaVecMaskALTIVEC< double > operator<(const InaVecALTIVEC< double >& inVec1, const InaVecALTIVEC< double >& inVec2) {
    return InaVecALTIVEC< double >::IsLowerMask(inVec1, inVec2);
}

inline InaVecMaskALTIVEC< double > operator<=(const InaVecALTIVEC< double >& inVec1, const InaVecALTIVEC< double >& inVec2) {
    return InaVecALTIVEC< double >::IsLowerOrEqualMask(inVec1, inVec2);
}

inline InaVecMaskALTIVEC< double > operator>(const InaVecALTIVEC< double >& inVec1, const InaVecALTIVEC< double >& inVec2) {
    return InaVecALTIVEC< double >::IsGreaterMask(inVec1, inVec2);
}

inline InaVecMaskALTIVEC< double > operator>=(const InaVecALTIVEC< double >& inVec1, const InaVecALTIVEC< double >& inVec2) {
    return InaVecALTIVEC< double >::IsGreaterOrEqualMask(inVec1, inVec2);
}

inline InaVecMaskALTIVEC< double > operator==(const InaVecALTIVEC< double >& inVec1, const InaVecALTIVEC< double >& inVec2) {
    return InaVecALTIVEC< double >::IsEqualMask(inVec1, inVec2);
}

inline InaVecMaskALTIVEC< double > operator!=(const InaVecALTIVEC< double >& inVec1, const InaVecALTIVEC< double >& inVec2) {
    return InaVecALTIVEC< double >::IsNotEqualMask(inVec1, inVec2);
}

#endif
