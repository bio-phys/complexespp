///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#include "InastempConfig.h"

#include "@TYPE@/InaVec@TYPE@Double.hpp"
#include "@TYPE@/InaVec@TYPE@Float.hpp"

#include "flops-test-all.hpp"
#include "FLOPS/InaVecFLOPS.hpp"

int main() {
    // clang-format off
    FlopsTestAll<InaVecFLOPS<InaVec@TYPE@<double> >> testerDoubleFlops;
    FlopsTestAll<InaVecFLOPS<InaVec@TYPE@<float> >> testerSingleFlops;
    // clang-format on
    return testerDoubleFlops.Run() + testerSingleFlops.Run();
}
