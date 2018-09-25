///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#include "InastempConfig.h"

#include "@TYPE@/InaVec@TYPE@Double.hpp"
#include "@TYPE@/InaVec@TYPE@Float.hpp"

#include "core-test-all.hpp"
#include "FLOPS/InaVecFLOPS.hpp"

int main() {
    // clang-format off
    TestAll<InaVec@TYPE@<double> > testerDouble;
    TestAll<InaVec@TYPE@<float> > testerSingle;

    TestAll<InaVecFLOPS<InaVec@TYPE@<double> >> testerDoubleFlops;
    TestAll<InaVecFLOPS<InaVec@TYPE@<float> >> testerSingleFlops;
    // clang-format on
    return testerDouble.Run() + testerSingle.Run() + testerDoubleFlops.Run() + testerSingleFlops.Run();
}
