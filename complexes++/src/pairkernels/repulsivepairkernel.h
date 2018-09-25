#ifndef REPULSIVEPAIRKERNEL_H
#define REPULSIVEPAIRKERNEL_H

#include <memory>

#include "abstractpairkernel.h"
#include "domains/abstractdomain.h"
#include "energy/forcefield.h"

#include "vectorization/energycomputervec.h"
#include "vectorization/repulsivekernelvec.h"

namespace pairkernels {

class RepulsivePairKernel : public AbstractPairKernel {
 public:
  static const int RepulsiveNbContributions = 1;

  static std::unique_ptr<AbstractPairKernel> BuildRepulsivePairKernel() {
    return std::make_unique<RepulsivePairKernel>();
  }

  explicit RepulsivePairKernel()
      : AbstractPairKernel(RepulsiveNbContributions) {}

  void compute(energy::EnergyMatrixBuffer<>::Accesser inView,
               const domains::AbstractDomain& inDom1,
               const domains::AbstractDomain& inDom2, const util::rvec& box,
               const energy::ForceField& forceField,
               const std::pair<int, int>& inInterval1,
               const std::pair<int, int>& inInterval2) const final {
    const double debyeLength(forceField.debyeLength());
    const double dielectricConstant(forceField.dielectricConstant());
    simd::EnergyComputerVecBest<simd::RepulsiveKernelVec, double> computer(
        debyeLength, dielectricConstant, forceField);
    computer.disctinctDomainsIntervals(inDom1, inDom2, box, inInterval1,
                                       inInterval2);
    // Contributions/energy are accumulated by the kernel retreive them and add
    // them to the view
    const double en = computer.template getContribution<
        simd::RepulsiveKernelVecContributions::En>();
    inView.addContribution(getOffsetContributions(), 0, en);
  }

  static std::string GetName() { return "Repulsive"; }

  std::string type() const override { return GetName(); }

  std::vector<std::string> getExtraContributionLabels() const final {
    std::vector<std::string> labels;
    labels.emplace_back("Repulsive");
    return labels;
  }
};
}  // namespace pairkernels

#endif
