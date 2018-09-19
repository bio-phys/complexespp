#ifndef LJHCUTOFFPAIRKERNEL_H
#define LJHCUTOFFPAIRKERNEL_H

#include <memory>
#include <vector>

#include "abstractpairkernel.h"
#include "domains/abstractdomain.h"
#include "energy/forcefield.h"

#include "vectorization/energycomputervec.h"
#include "vectorization/ljhcutoffkernelvec.h"

namespace pairkernels {

class LJHCutoffPairKernel : public AbstractPairKernel {
 public:
  static const int LJHNbContributions = 2;

  static std::unique_ptr<AbstractPairKernel> BuildLJHCutoffPairKernel() {
    return std::make_unique<LJHCutoffPairKernel>();
  }

  explicit LJHCutoffPairKernel() : AbstractPairKernel(LJHNbContributions) {}

  void compute(energy::EnergyMatrixBuffer<>::Accesser inView,
               const domains::AbstractDomain& inDom1,
               const domains::AbstractDomain& inDom2, const util::rvec& box,
               const energy::ForceField& forceField,
               const std::pair<int, int>& inInterval1,
               const std::pair<int, int>& inInterval2) const final {
    const double debyeLength(forceField.debyeLength());
    const double dielectricConstant(forceField.dielectricConstant());
    simd::EnergyComputerVecBest<simd::LJHCutoffKernelVec, double> computer(
        debyeLength, dielectricConstant, forceField);
    computer.disctinctDomainsIntervals(inDom1, inDom2, box, inInterval1,
                                       inInterval2);
    // Contributions/energy are accumulated by the kernel
    // retreive them and add them to the view
    const double lennardJones =
        computer.getContribution<simd::LJHCutoffKernelVecContributions::LJH>();
    const double debyeHueckel =
        computer.getContribution<simd::LJHCutoffKernelVecContributions::DH>();

    inView.addContribution(getOffsetContributions(), 0, lennardJones);
    inView.addContribution(getOffsetContributions(), 1, debyeHueckel);
  }

  static std::string GetName() { return "LJHCutoff"; }

  std::string type() const override { return GetName(); }

  std::vector<std::string> getExtraContributionLabels() const final {
    std::vector<std::string> labels;
    labels.emplace_back("LJHCutoff-lennardJones");
    labels.emplace_back("LJHCutoff-debyeHueckel");
    return labels;
  }
};
}  // namespace pairkernels

#endif
