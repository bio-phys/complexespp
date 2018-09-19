///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////

#include "InastempConfig.h"

#include "@TYPE@/InaVec@TYPE@Double.hpp"
#include "@TYPE@/InaVec@TYPE@Float.hpp"

#include "UTester.hpp"

#include <cmath>
#include <cstring>

template < class VecType >
class TestAll : public UTester< TestAll< VecType > > {
    using Parent = UTester< TestAll< VecType > >;

    using RealType = typename VecType::RealType;

    void TestBasic() {
        RealType inputVal[VecType::VecLength];
        for (int idx = 0; idx < VecType::VecLength; ++idx) {
            inputVal[idx] = RealType(idx);
        }
        const VecType inputVec(inputVal);

        {
            const VecType res0 = VecType::If(inputVec.isZeroMask()).Then([&]() {
                return RealType(1);
            });

            RealType outputVal[VecType::VecLength];
            res0.storeInArray(outputVal);

            for (int idx = 0; idx < VecType::VecLength; ++idx) {
                if (inputVal[idx] == 0) {
                    UASSERTEEQUAL(outputVal[idx], RealType(1));
                } else {
                    UASSERTEEQUAL(outputVal[idx], RealType(0));
                }
            }
        }
        {
            const VecType res0 = VecType::If(inputVec.isZeroMask())
                                     .Then([&]() {
                                         return RealType(-1);
                                     })
                                     .Else([&]() {
                                         return RealType(99);
                                     });

            RealType outputVal[VecType::VecLength];
            res0.storeInArray(outputVal);

            for (int idx = 0; idx < VecType::VecLength; ++idx) {
                if (inputVal[idx] == 0) {
                    UASSERTEEQUAL(outputVal[idx], RealType(-1));
                } else {
                    UASSERTEEQUAL(outputVal[idx], RealType(99));
                }
            }
        }
        {
            const VecType res0 = VecType::If(inputVec.isZeroMask())
                                     .Then([&]() {
                                         return RealType(1);
                                     })
                                     .ElseIf(VecType::IsLowerOrEqualMask(1, inputVec))
                                     .Then([&]() {
                                         return RealType(2);
                                     })
                                     .Else([&]() {
                                         return RealType(3);
                                     });

            RealType outputVal[VecType::VecLength];
            res0.storeInArray(outputVal);

            for (int idx = 0; idx < VecType::VecLength; ++idx) {
                if (inputVal[idx] == 0) {
                    UASSERTEEQUAL(outputVal[idx], RealType(1));
                } else if (1 <= inputVal[idx]) {
                    UASSERTEEQUAL(outputVal[idx], RealType(2));
                } else {
                    UASSERTEEQUAL(outputVal[idx], RealType(3));
                }
            }
        }
        {
            const VecType res0 = VecType::If(inputVec.isNotZeroMask())
                                     .Then([&]() {
                                         return RealType(1);
                                     })
                                     .ElseIf(VecType::IsGreaterOrEqualMask(1, inputVec))
                                     .Then([&]() {
                                         return RealType(2);
                                     })
                                     .Else([&]() {
                                         return RealType(3);
                                     });

            RealType outputVal[VecType::VecLength];
            res0.storeInArray(outputVal);

            for (int idx = 0; idx < VecType::VecLength; ++idx) {
                if (inputVal[idx] != 0) {
                    UASSERTEEQUAL(outputVal[idx], RealType(1));
                } else if (1 >= inputVal[idx]) {
                    UASSERTEEQUAL(outputVal[idx], RealType(2));
                } else {
                    UASSERTEEQUAL(outputVal[idx], RealType(3));
                }
            }
        }
        {
            const VecType res0 = VecType::If(VecType::IsEqualMask(100,
                                                                inputVec))
                                     .Then([&]() {
                                         return RealType(1);
                                     })
                                     .ElseIf(VecType::IsGreaterOrEqualMask(1, inputVec))
                                     .Then([&]() {
                                         return RealType(2);
                                     })
                                     .Else([&]() {
                                         return RealType(3);
                                     });

            RealType outputVal[VecType::VecLength];
            res0.storeInArray(outputVal);

            for (int idx = 0; idx < VecType::VecLength; ++idx) {
                if (inputVal[idx] == 100) {
                    UASSERTEEQUAL(outputVal[idx], RealType(1));
                } else if (1 >= inputVal[idx]) {
                    UASSERTEEQUAL(outputVal[idx], RealType(2));
                } else {
                    UASSERTEEQUAL(outputVal[idx], RealType(3));
                }
            }
        }

        {
            const VecType res0 = VecType::IfElse(inputVec.isZeroMask(),
                                                   -1,
                                                   99);

            RealType outputVal[VecType::VecLength];
            res0.storeInArray(outputVal);

            for (int idx = 0; idx < VecType::VecLength; ++idx) {
                if (inputVal[idx] == 0) {
                    UASSERTEEQUAL(outputVal[idx], RealType(-1));
                } else {
                    UASSERTEEQUAL(outputVal[idx], RealType(99));
                }
            }
        }
        {
            const VecType res0 = VecType::IfElse(inputVec.isNotZeroMask(),
                                                   -1,
                                                   99);

            RealType outputVal[VecType::VecLength];
            res0.storeInArray(outputVal);

            for (int idx = 0; idx < VecType::VecLength; ++idx) {
                if (inputVal[idx] != 0) {
                    UASSERTEEQUAL(outputVal[idx], RealType(-1));
                } else {
                    UASSERTEEQUAL(outputVal[idx], RealType(99));
                }
            }
        }

        {
            const VecType res0 = VecType::IfTrue(inputVec.isZeroMask(),
                                                   -1);

            RealType outputVal[VecType::VecLength];
            res0.storeInArray(outputVal);

            for (int idx = 0; idx < VecType::VecLength; ++idx) {
                if (inputVal[idx] == 0) {
                    UASSERTEEQUAL(outputVal[idx], RealType(-1));
                } else {
                    UASSERTEEQUAL(outputVal[idx], RealType(0));
                }
            }
        }
        {
            const VecType res0 = VecType::IfTrue(inputVec.isNotZeroMask(),
                                                   -1);

            RealType outputVal[VecType::VecLength];
            res0.storeInArray(outputVal);

            for (int idx = 0; idx < VecType::VecLength; ++idx) {
                if (inputVal[idx] != 0) {
                    UASSERTEEQUAL(outputVal[idx], RealType(-1));
                } else {
                    UASSERTEEQUAL(outputVal[idx], RealType(0));
                }
            }
        }
    }


    void TestOperator() {
        RealType inputVal[VecType::VecLength];
        for (int idx = 0; idx < VecType::VecLength; ++idx) {
            inputVal[idx] = RealType(idx);
        }
        const VecType inputVec(inputVal);

        {
            const VecType res0 = VecType::If(inputVec == 0).Then([&]() {
                return RealType(1);
            });

            RealType outputVal[VecType::VecLength];
            res0.storeInArray(outputVal);

            for (int idx = 0; idx < VecType::VecLength; ++idx) {
                if (inputVal[idx] == 0) {
                    UASSERTEEQUAL(outputVal[idx], RealType(1));
                } else {
                    UASSERTEEQUAL(outputVal[idx], RealType(0));
                }
            }
        }
        {
            const VecType res0 = VecType::If(inputVec == 0)
                                     .Then([&]() {
                                         return RealType(-1);
                                     })
                                     .Else([&]() {
                                         return RealType(99);
                                     });

            RealType outputVal[VecType::VecLength];
            res0.storeInArray(outputVal);

            for (int idx = 0; idx < VecType::VecLength; ++idx) {
                if (inputVal[idx] == 0) {
                    UASSERTEEQUAL(outputVal[idx], RealType(-1));
                } else {
                    UASSERTEEQUAL(outputVal[idx], RealType(99));
                }
            }
        }
        {
            const VecType res0 = VecType::If(inputVec == 0)
                                     .Then([&]() {
                                         return RealType(1);
                                     })
                                     .ElseIf(inputVec >= 1)
                                     .Then([&]() {
                                         return RealType(2);
                                     })
                                     .Else([&]() {
                                         return RealType(3);
                                     });

            RealType outputVal[VecType::VecLength];
            res0.storeInArray(outputVal);

            for (int idx = 0; idx < VecType::VecLength; ++idx) {
                if (inputVal[idx] == 0) {
                    UASSERTEEQUAL(outputVal[idx], RealType(1));
                } else if (1 <= inputVal[idx]) {
                    UASSERTEEQUAL(outputVal[idx], RealType(2));
                } else {
                    UASSERTEEQUAL(outputVal[idx], RealType(3));
                }
            }
        }
        {
            const VecType res0 = VecType::If(inputVec != 0)
                                     .Then([&]() {
                                         return RealType(1);
                                     })
                                     .ElseIf(inputVec <= 1)
                                     .Then([&]() {
                                         return RealType(2);
                                     })
                                     .Else([&]() {
                                         return RealType(3);
                                     });

            RealType outputVal[VecType::VecLength];
            res0.storeInArray(outputVal);

            for (int idx = 0; idx < VecType::VecLength; ++idx) {
                if (inputVal[idx] != 0) {
                    UASSERTEEQUAL(outputVal[idx], RealType(1));
                } else if (1 >= inputVal[idx]) {
                    UASSERTEEQUAL(outputVal[idx], RealType(2));
                } else {
                    UASSERTEEQUAL(outputVal[idx], RealType(3));
                }
            }
        }
        {
            const VecType res0 = VecType::If(inputVec == 100)
                                     .Then([&]() {
                                         return RealType(1);
                                     })
                                     .ElseIf(inputVec <= 1)
                                     .Then([&]() {
                                         return RealType(2);
                                     })
                                     .Else([&]() {
                                         return RealType(3);
                                     });

            RealType outputVal[VecType::VecLength];
            res0.storeInArray(outputVal);

            for (int idx = 0; idx < VecType::VecLength; ++idx) {
                if (inputVal[idx] == 100) {
                    UASSERTEEQUAL(outputVal[idx], RealType(1));
                } else if (1 >= inputVal[idx]) {
                    UASSERTEEQUAL(outputVal[idx], RealType(2));
                } else {
                    UASSERTEEQUAL(outputVal[idx], RealType(3));
                }
            }
        }

        {
            const VecType res0 = VecType::IfElse(inputVec == 0, -1, 99);

            RealType outputVal[VecType::VecLength];
            res0.storeInArray(outputVal);

            for (int idx = 0; idx < VecType::VecLength; ++idx) {
                if (inputVal[idx] == 0) {
                    UASSERTEEQUAL(outputVal[idx], RealType(-1));
                } else {
                    UASSERTEEQUAL(outputVal[idx], RealType(99));
                }
            }
        }
        {
            const VecType res0 = VecType::IfElse(inputVec != 0, -1, 99);

            RealType outputVal[VecType::VecLength];
            res0.storeInArray(outputVal);

            for (int idx = 0; idx < VecType::VecLength; ++idx) {
                if (inputVal[idx] != 0) {
                    UASSERTEEQUAL(outputVal[idx], RealType(-1));
                } else {
                    UASSERTEEQUAL(outputVal[idx], RealType(99));
                }
            }
        }

        {
            const VecType res0 = VecType::IfTrue(inputVec == 0, -1);

            RealType outputVal[VecType::VecLength];
            res0.storeInArray(outputVal);

            for (int idx = 0; idx < VecType::VecLength; ++idx) {
                if (inputVal[idx] == 0) {
                    UASSERTEEQUAL(outputVal[idx], RealType(-1));
                } else {
                    UASSERTEEQUAL(outputVal[idx], RealType(0));
                }
            }
        }
        {
            const VecType res0 = VecType::IfTrue(inputVec != 0, -1);

            RealType outputVal[VecType::VecLength];
            res0.storeInArray(outputVal);

            for (int idx = 0; idx < VecType::VecLength; ++idx) {
                if (inputVal[idx] != 0) {
                    UASSERTEEQUAL(outputVal[idx], RealType(-1));
                } else {
                    UASSERTEEQUAL(outputVal[idx], RealType(0));
                }
            }
        }
    }


    void TestOperatorNoRet() {
        RealType inputVal[VecType::VecLength];
        for (int idx = 0; idx < VecType::VecLength; ++idx) {
            inputVal[idx] = RealType(idx);
        }
        const VecType inputVec(inputVal);

        {
            const VecType res0 = VecType::If(inputVec == 0).Then(1);

            RealType outputVal[VecType::VecLength];
            res0.storeInArray(outputVal);

            for (int idx = 0; idx < VecType::VecLength; ++idx) {
                if (inputVal[idx] == 0) {
                    UASSERTEEQUAL(outputVal[idx], RealType(1));
                } else {
                    UASSERTEEQUAL(outputVal[idx], RealType(0));
                }
            }
        }
        {
            const VecType res0 = VecType::If(inputVec == 0).Then(-1).Else(99);

            RealType outputVal[VecType::VecLength];
            res0.storeInArray(outputVal);

            for (int idx = 0; idx < VecType::VecLength; ++idx) {
                if (inputVal[idx] == 0) {
                    UASSERTEEQUAL(outputVal[idx], RealType(-1));
                } else {
                    UASSERTEEQUAL(outputVal[idx], RealType(99));
                }
            }
        }
        {
            const VecType res0 = VecType::If(inputVec == 0).Then(1).ElseIf(inputVec >= 1).Then(2).Else(3);

            RealType outputVal[VecType::VecLength];
            res0.storeInArray(outputVal);

            for (int idx = 0; idx < VecType::VecLength; ++idx) {
                if (inputVal[idx] == 0) {
                    UASSERTEEQUAL(outputVal[idx], RealType(1));
                } else if (1 <= inputVal[idx]) {
                    UASSERTEEQUAL(outputVal[idx], RealType(2));
                } else {
                    UASSERTEEQUAL(outputVal[idx], RealType(3));
                }
            }
        }
        {
            const VecType res0 = VecType::If(inputVec != 0).Then(1).ElseIf(inputVec <= 1).Then(2).Else(3);

            RealType outputVal[VecType::VecLength];
            res0.storeInArray(outputVal);

            for (int idx = 0; idx < VecType::VecLength; ++idx) {
                if (inputVal[idx] != 0) {
                    UASSERTEEQUAL(outputVal[idx], RealType(1));
                } else if (1 >= inputVal[idx]) {
                    UASSERTEEQUAL(outputVal[idx], RealType(2));
                } else {
                    UASSERTEEQUAL(outputVal[idx], RealType(3));
                }
            }
        }
        {
            const VecType res0 = VecType::If(inputVec == 100).Then(1).ElseIf(inputVec <= 1).Then(2).Else(3);

            RealType outputVal[VecType::VecLength];
            res0.storeInArray(outputVal);

            for (int idx = 0; idx < VecType::VecLength; ++idx) {
                if (inputVal[idx] == 100) {
                    UASSERTEEQUAL(outputVal[idx], RealType(1));
                } else if (1 >= inputVal[idx]) {
                    UASSERTEEQUAL(outputVal[idx], RealType(2));
                } else {
                    UASSERTEEQUAL(outputVal[idx], RealType(3));
                }
            }
        }
        {
            UASSERTETRUE((inputVec == 100.) == (inputVec == 100.));
            UASSERTETRUE((inputVec != 100.) != (inputVec == 100.));
            UASSERTETRUE((inputVec != 100.) == (inputVec != 100.));
            UASSERTETRUE((VecType(0.) == VecType(1.)).isAllTrue() == false);
            UASSERTETRUE((VecType(0.) == VecType(1.)).isAllFalse() == true);
            UASSERTETRUE((VecType(0.) == VecType(0.)).isAllFalse() == false);
            UASSERTETRUE((VecType(0.) == VecType(0.)).isAllTrue() == true);

            RealType values10[VecType::VecLength] = {0};
            values10[0] = 1;
            VecType vecValues10(values10);

            UASSERTETRUE((vecValues10 == 100.) == (vecValues10 == 100.));
            UASSERTETRUE((vecValues10 != 100.) != (vecValues10 == 100.));
            UASSERTETRUE((vecValues10 != 100.) == (vecValues10 != 100.));
            UASSERTETRUE((vecValues10 == VecType(1.)).isAllTrue() == false || VecType::VecLength == 1);
            UASSERTETRUE((vecValues10 == VecType(1.)).isAllFalse() == false);
            UASSERTETRUE((vecValues10 == VecType(0.)).isAllFalse() == false || VecType::VecLength == 1);
            UASSERTETRUE((vecValues10 == VecType(0.)).isAllTrue() == false);

            UASSERTETRUE((vecValues10 == vecValues10).isAllTrue() == true);
            UASSERTETRUE((vecValues10 != vecValues10).isAllFalse() == true);
        }
    }

    void SetTests() {
        Parent::AddTest(&TestAll::TestBasic, "Basic ifelse tests for vec type");
        Parent::AddTest(&TestAll::TestOperator, "Basic ifelse tests for vec type but with overloaded operators");
        Parent::AddTest(&TestAll::TestOperatorNoRet, "Basic ifelse tests for vec type but with overloaded operators"
                                                     "and no return in statement");
    }
};


int main() {
    // clang-format off
    TestAll< InaVec@TYPE@<double> > testerDouble;
    TestAll< InaVec@TYPE@<float> > testerSingle;
    // clang-format on
    return testerDouble.Run() + testerSingle.Run();
}
