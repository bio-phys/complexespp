///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////

#include "InastempConfig.h"


#include "@TYPE@/InaVec@TYPE@Double.hpp"
#include "@TYPE@/InaVec@TYPE@Float.hpp"


int main() {
    // clang-format off
    {
        InaVec@TYPE@<double> floatVec;
        InaVec@TYPE@<double> vecOne = 1;
        typename InaVec@TYPE@<double>::MaskType msk = (1 < vecOne);
        msk = msk;
    }
    {
        InaVec@TYPE@<float> floatVec;
        InaVec@TYPE@<float> vecOne = 1;
        typename InaVec@TYPE@<float>::MaskType msk = (1 < vecOne);
        msk = msk;
    }
    // clang-format on

    return 0;
}
