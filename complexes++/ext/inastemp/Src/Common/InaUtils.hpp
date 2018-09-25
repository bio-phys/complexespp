///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef UTILS_HPP
#define UTILS_HPP

// for std::size_t
#include <cstdio>

namespace InaUtils {

template < class TypeSrc, class TypeDest >
inline TypeDest ConvertBits(TypeSrc src) {
    union InaBits {
        TypeSrc src;
        TypeDest dest;
    };
    InaBits bits;
    bits.src = src;
    return bits.dest;
}

template < class VecType >
inline VecType FastPow(VecType base, std::size_t power) {
    VecType res = VecType(1.);

    while(power) {
        if(1 & power) {
            res *= base;
        }
        base *= base;
        power >>= 1;
    }

    return res;
}

inline std::size_t FastPowNbMul(std::size_t power) {
    std::size_t nbMul = 0;

    while(power) {
        if(1 & power) {
            nbMul += 1;
        }
        nbMul += 1;
        power >>= 1;
    }

    return nbMul;
}
}

#endif // UTILS_HPP
