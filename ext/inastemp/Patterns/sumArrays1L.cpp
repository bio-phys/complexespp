///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#include "InastempConfig.h"
#include "SCALAR/InaVecSCALARDouble.hpp"

#include <cassert>
#include <iostream>
#include <memory>

// Here we illustrate how to compute one loop based kernel.
// We compute res[i] = (a[i] + b[i]) for all i

template < class VecType >
void SumArrays(double* __restrict__ dest, const double* __restrict__ src1,
               const double* __restrict__ src2, const size_t nbToProceed) {

    for (size_t idx = 0; idx < nbToProceed; idx += VecType::VecLength) {
        const VecType v1(&src1[idx]);
        const VecType v2(&src2[idx]);
        const VecType res = v1 + v2;
        res.storeInArray(&dest[idx]);
        // In one line:
        //VecType(VecType(&src1[idx]) + VecType(&src2[idx]).storeInArray(&dest[idx]);
    }
}

void Pattern1D(double* __restrict__ dest, const double* __restrict__ src1,
               const double* __restrict__ src2, const size_t nbToProceed) {
    const size_t nbItemsVectorized = (nbToProceed / InaVecBestType<double>::VecLength) * InaVecBestType<double>::VecLength;
    // First with the best vectorizer
    SumArrays< InaVecBestType<double> >(dest, src1, src2, nbItemsVectorized);
    // Do the rest with scalar
    SumArrays< InaVecSCALAR<double> >(&dest[nbItemsVectorized], &src1[nbItemsVectorized],
                                    &src2[nbItemsVectorized], nbToProceed - nbItemsVectorized);
}

int main() {
    std::cout << "The best vectorizer computes " << InaVecBestType<double>::VecLength << " double values together." << std::endl;

    // Test size
    const size_t Size = 10000;

    // Allocate three arrays for testing
    std::unique_ptr< double[] > arraySrc1(new double[Size]);
    std::unique_ptr< double[] > arraySrc2(new double[Size]);
    std::unique_ptr< double[] > arrayRes(new double[Size]);

    // Fill source arrays with dumb values
    for (size_t idx = 0; idx < Size; ++idx) {
        arraySrc1[idx] = static_cast<double>(idx);
        arraySrc2[idx] = static_cast<double>(idx);
        arrayRes[idx]  = 0;
    }

    // Call kernel
    Pattern1D(arrayRes.get(), arraySrc1.get(), arraySrc2.get(), Size);

    // Just to check results (we compare bits here)
    for (size_t idx = 0; idx < Size; ++idx) {
        assert(arrayRes[idx] == static_cast<double>(2 * idx));
    }

    std::cout << "Done." << std::endl;

    return 0;
}
