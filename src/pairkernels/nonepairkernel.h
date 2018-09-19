#ifndef NONEPAIRKERNEL_H
#define NONEPAIRKERNEL_H

#include <memory>

#include "abstractpairkernel.h"
#include "domains/abstractdomain.h"
#include "energy/forcefield.h"

namespace pairkernels {

class NonePairKernel : public AbstractPairKernel {
 public:
  using AbstractPairKernel::AbstractPairKernel;

  static std::unique_ptr<AbstractPairKernel> BuildNonePairKernel() {
    return std::make_unique<NonePairKernel>(0);
  }

  void compute(energy::EnergyMatrixBuffer<>::Accesser /*inView*/,
               const domains::AbstractDomain& /*inDom1*/,
               const domains::AbstractDomain& /*inDom2*/,
               const util::rvec& /*box*/,
               const energy::ForceField& /*forcefield*/,
               const std::pair<int, int>& /*inInterval1*/,
               const std::pair<int, int>& /*inInterval2*/) const final {}

  static std::string GetName() { return "None"; }

  std::string type() const override { return GetName(); }
};
}  // namespace pairkernels

#endif
