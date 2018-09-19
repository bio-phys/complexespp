///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef INAIFELSE_HPP
#define INAIFELSE_HPP

#include <functional>
#include <type_traits>

template < class VecType >
class InaIfElse {
public:
    using MaskType = typename VecType::MaskType;

    class ThenClass;

    class ElseIfClass {
    protected:
        MaskType cumulatedTest;
        VecType currentValue;

    public:
        inline ElseIfClass(const MaskType& inCumulatedTest,
                    const VecType& inCurrentValue)
        : cumulatedTest(inCumulatedTest), currentValue(inCurrentValue) {
        }

        inline ThenClass ElseIf(const MaskType& inTestRes) {
            return ThenClass(inTestRes, cumulatedTest, currentValue);
        }

        template<class FunctionType>
        inline typename std::enable_if<!std::is_convertible<FunctionType,VecType>::value, VecType>::type
             Else(FunctionType snippet) {
            const VecType snippetRes = static_cast<VecType>(snippet());
            return Else(snippetRes);
        }

        inline VecType Else(const VecType& rawValue) {
            currentValue = VecType::BitsOr(VecType::IfFalse(cumulatedTest, rawValue), currentValue);
            return currentValue;
        }

        inline operator VecType() const {
            return currentValue;
        }
    };

    class ThenClass {
    protected:
        MaskType onGoingTest;
        MaskType cumulatedTest;
        VecType currentValue;

    public:
        inline ThenClass(const MaskType& inOnGoingTest,
                  const MaskType& inCumulatedTest, const VecType& inCurrentValue)
        : onGoingTest(inOnGoingTest),
          cumulatedTest(inCumulatedTest), currentValue(inCurrentValue) {
        }

        template<class FunctionType>
        inline typename std::enable_if<!std::is_convertible<FunctionType,VecType>::value, ElseIfClass>::type
                Then(FunctionType snippet) {
            const VecType snippetRes = static_cast<VecType>(snippet());
            return Then(snippetRes);
        }

        inline ElseIfClass Then(const VecType& rawValue) {
            currentValue = VecType::BitsOr(VecType::IfFalse(cumulatedTest,
                                                   VecType::IfTrue(onGoingTest, rawValue)),
                                             currentValue);
            cumulatedTest = MaskType::Or(cumulatedTest, onGoingTest);
            return ElseIfClass(cumulatedTest, currentValue);
        }

        inline operator VecType() const {
            return currentValue;
        }
    };

    class IfClass {
    public:
        inline IfClass() {
        }

        inline ThenClass If(const MaskType& inTestRes) {
            return ThenClass(inTestRes, MaskType(false), VecType::GetZero());
        }
    };
};

#endif
