///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef FLOPSTESTALL_HPP
#define FLOPSTESTALL_HPP

#include "InastempConfig.h"
#include "UTester.hpp"

#include <cmath>
#include <cstring>

template < class VecType >
class FlopsTestAll : public UTester< FlopsTestAll< VecType > > {
    using Parent = UTester< FlopsTestAll< VecType > >;

    using RealType = typename VecType::RealType;
    using MaskType = typename VecType::MaskType;

    void TestBasic() {
        UASSERTEEQUAL(VecType::GetFlopsStats().getMulOp(), VecType::VecLength * size_t(0));
        UASSERTEEQUAL(VecType::GetFlopsStats().getDivOp(), VecType::VecLength * size_t(0));
        UASSERTEEQUAL(VecType::GetFlopsStats().getAddOp(), VecType::VecLength * size_t(0));
        UASSERTEEQUAL(VecType::GetFlopsStats().getSubOp(), VecType::VecLength * size_t(0));
        UASSERTEEQUAL(VecType::GetFlopsStats().getRsqrt(), VecType::VecLength * size_t(0));
        UASSERTEEQUAL(VecType::GetFlopsStats().getSqrt(), VecType::VecLength * size_t(0));


        VecType a = 1;
        {
            VecType res = a + a;
            UASSERTEEQUAL(VecType::GetFlopsStats().getMulOp(), VecType::VecLength * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getDivOp(), VecType::VecLength * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getAddOp(), VecType::VecLength * size_t(1));
            UASSERTEEQUAL(VecType::GetFlopsStats().getSubOp(), VecType::VecLength * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getRsqrt(), VecType::VecLength * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getSqrt(), VecType::VecLength * size_t(0));

            res += a;
            UASSERTEEQUAL(VecType::GetFlopsStats().getMulOp(), VecType::VecLength * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getDivOp(), VecType::VecLength * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getAddOp(), VecType::VecLength * size_t(2));
            UASSERTEEQUAL(VecType::GetFlopsStats().getSubOp(), VecType::VecLength * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getRsqrt(), VecType::VecLength * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getSqrt(), VecType::VecLength * size_t(0));
        }

        VecType::ResetFlopsStats();
        {
            VecType res = a * a;
            UASSERTEEQUAL(VecType::GetFlopsStats().getMulOp(), VecType::VecLength * size_t(1));
            UASSERTEEQUAL(VecType::GetFlopsStats().getDivOp(), VecType::VecLength * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getAddOp(), VecType::VecLength * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getSubOp(), VecType::VecLength * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getRsqrt(), VecType::VecLength * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getSqrt(), VecType::VecLength * size_t(0));

            res *= a;
            UASSERTEEQUAL(VecType::GetFlopsStats().getMulOp(), VecType::VecLength * size_t(2));
            UASSERTEEQUAL(VecType::GetFlopsStats().getDivOp(), VecType::VecLength * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getAddOp(), VecType::VecLength * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getSubOp(), VecType::VecLength * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getRsqrt(), VecType::VecLength * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getSqrt(), VecType::VecLength * size_t(0));
        }

        VecType::ResetFlopsStats();
        {
            VecType res = a / a;
            UASSERTEEQUAL(VecType::GetFlopsStats().getMulOp(), VecType::VecLength * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getDivOp(), VecType::VecLength * size_t(1));
            UASSERTEEQUAL(VecType::GetFlopsStats().getAddOp(), VecType::VecLength * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getSubOp(), VecType::VecLength * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getRsqrt(), VecType::VecLength * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getSqrt(), VecType::VecLength * size_t(0));

            res /= a;
            UASSERTEEQUAL(VecType::GetFlopsStats().getMulOp(), VecType::VecLength * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getDivOp(), VecType::VecLength * size_t(2));
            UASSERTEEQUAL(VecType::GetFlopsStats().getAddOp(), VecType::VecLength * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getSubOp(), VecType::VecLength * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getRsqrt(), VecType::VecLength * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getSqrt(), VecType::VecLength * size_t(0));
        }

        VecType::ResetFlopsStats();
        {
            VecType res = a - a;
            UASSERTEEQUAL(VecType::GetFlopsStats().getMulOp(), VecType::VecLength * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getDivOp(), VecType::VecLength * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getAddOp(), VecType::VecLength * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getSubOp(), VecType::VecLength * size_t(1));
            UASSERTEEQUAL(VecType::GetFlopsStats().getRsqrt(), VecType::VecLength * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getSqrt(), VecType::VecLength * size_t(0));

            res -= a;
            UASSERTEEQUAL(VecType::GetFlopsStats().getMulOp(), VecType::VecLength * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getDivOp(), VecType::VecLength * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getAddOp(), VecType::VecLength * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getSubOp(), VecType::VecLength * size_t(2));
            UASSERTEEQUAL(VecType::GetFlopsStats().getRsqrt(), VecType::VecLength * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getSqrt(), VecType::VecLength * size_t(0));
        }


        VecType::ResetFlopsStats();
        {
            VecType res = (a * a) + (a / a) - (a + a) * (a - a) / a;

            UASSERTEEQUAL(VecType::GetFlopsStats().getMulOp(), VecType::VecLength * size_t(2));
            UASSERTEEQUAL(VecType::GetFlopsStats().getDivOp(), VecType::VecLength * size_t(2));
            UASSERTEEQUAL(VecType::GetFlopsStats().getAddOp(), VecType::VecLength * size_t(2));
            UASSERTEEQUAL(VecType::GetFlopsStats().getSubOp(), VecType::VecLength * size_t(2));
            UASSERTEEQUAL(VecType::GetFlopsStats().getRsqrt(), VecType::VecLength * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getSqrt(), VecType::VecLength * size_t(0));

            res = VecType(0.);

            UASSERTEEQUAL(VecType::GetFlopsStats().getMulOp(), VecType::VecLength * size_t(2));
            UASSERTEEQUAL(VecType::GetFlopsStats().getDivOp(), VecType::VecLength * size_t(2));
            UASSERTEEQUAL(VecType::GetFlopsStats().getAddOp(), VecType::VecLength * size_t(2));
            UASSERTEEQUAL(VecType::GetFlopsStats().getSubOp(), VecType::VecLength * size_t(2));
            UASSERTEEQUAL(VecType::GetFlopsStats().getRsqrt(), VecType::VecLength * size_t(0));
            UASSERTEEQUAL(VecType::GetFlopsStats().getSqrt(), VecType::VecLength * size_t(0));
        }
    }

    void SetTests() {
        Parent::AddTest(&FlopsTestAll::TestBasic, "Basic test for vec type");
    }
};

#endif
