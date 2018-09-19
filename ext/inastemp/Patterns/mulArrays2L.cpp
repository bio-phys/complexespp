///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#include "InastempConfig.h"
#include "SCALAR/InaVecSCALARDouble.hpp"

#include <cassert>
#include <iostream>
#include <memory>

// Here we illustrate how to compute pair-wise interactions.
// We compute Sum(a[i] * b[j]) for all i,j possible


template < class VecType >
double MulArrays(const double* __restrict__ src1, const size_t nbToProceed1,
                 const double* __restrict__ src2, const size_t nbToProceed2) {

    VecType sum = VecType::GetZero();

    for (size_t idx1 = 0; idx1 < nbToProceed1; ++idx1) {
        const VecType v1(src1[idx1]);

        for (size_t idx2 = 0; idx2 < nbToProceed2; idx2 += VecType::VecLength) {
            const VecType v2(&src2[idx2]);
            sum += v1 * v2;
            // In one line:
            // sum += v1 * VecType(&src2[idx]);
        }
    }
    return sum.horizontalSum();
}

double Pattern2D(const double* __restrict__ src1, const size_t nbToProceed1,
                 const double* __restrict__ src2, const size_t nbToProceed2) {
    const size_t nbItemsVectorized2 = (nbToProceed2 / InaVecBestType<double>::VecLength) * InaVecBestType<double>::VecLength;
    // First with the best vectorizer
    double res = MulArrays< InaVecBestType<double> >(src1, nbToProceed1, src2, nbItemsVectorized2);
    // Do the rest with scalar
    res += MulArrays< InaVecSCALAR<double> >(src1, nbToProceed1,
                                           &src2[nbItemsVectorized2], nbToProceed2 - nbItemsVectorized2);
    return res;
}

int main() {
    std::cout << "The best vectorizer computes " << InaVecBestType<double>::VecLength << " double values together." << std::endl;

    // Test size
    const size_t Size = 10000;

    // Allocate three arrays for testing
    std::unique_ptr< double[] > arraySrc1(new double[Size]);
    std::unique_ptr< double[] > arraySrc2(new double[Size]);

    // Fill source arrays with dumb values
    for (size_t idx = 0; idx < Size; ++idx) {
        arraySrc1[idx] = static_cast<double>(idx);
        arraySrc2[idx] = static_cast<double>(idx);
    }

    // Call kernel
    /*const double totalRes = */ Pattern2D(arraySrc1.get(), Size, arraySrc2.get(), Size);

    std::cout << "Done." << std::endl;

    return 0;
}
