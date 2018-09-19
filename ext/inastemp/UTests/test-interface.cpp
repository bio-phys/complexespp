///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////

#include "InastempConfig.h"

#include "Common/InaVecInterface.hpp"

#include "@TYPE@/InaVec@TYPE@Double.hpp"
#include "@TYPE@/InaVec@TYPE@Float.hpp"

#include "core-test-all.hpp"

int main() {
    // clang-format off
    {
        InaVecMaskInterface<InaVec@TYPE@<double>::MaskType> mask;
        InaVecInterface<InaVec@TYPE@<double>> vec;
    }
    {
        InaVecMaskInterface<InaVec@TYPE@<float>::MaskType> mask;
        InaVecInterface<InaVec@TYPE@<float>> vec;
    }
    TestAll<InaVecInterface<InaVec@TYPE@<double>>> testerDouble;
    TestAll<InaVecInterface<InaVec@TYPE@<float>>> testerSingle;
    // clang-format on
    return testerDouble.Run() + testerSingle.Run();
}
