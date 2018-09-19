#ifndef LJHPAIRKERNEL_H
#define LJHPAIRKERNEL_H

#include <memory>
#include <vector>

#include "abstractpairkernel.h"
#include "domains/abstractdomain.h"
#include "energy/forcefield.h"

#include "vectorization/energycomputervec.h"
#include "vectorization/ljhkernelvec.h"

namespace pairkernels {

class LJHPairKernel : public AbstractPairKernel {
 public:
  static const int LJHNbContributions = 2;

  static std::unique_ptr<AbstractPairKernel> BuildLJHPairKernel() {
    return std::make_unique<LJHPairKernel>();
  }

  explicit LJHPairKernel() : AbstractPairKernel(LJHNbContributions) {}

  void compute(energy::EnergyMatrixBuffer<>::Accesser inView,
               const domains::AbstractDomain& inDom1,
               const domains::AbstractDomain& inDom2, const util::rvec& box,
               const energy::ForceField& forceField,
               const std::pair<int, int>& inInterval1,
               const std::pair<int, int>& inInterval2) const final {
    const double debyeLength(forceField.debyeLength());
    const double dielectricConstant(forceField.dielectricConstant());
    simd::EnergyComputerVecBest<simd::LJHKernelVec, double> computer(
        debyeLength, dielectricConstant, forceField);
    computer.disctinctDomainsIntervals(inDom1, inDom2, box, inInterval1,
                                       inInterval2);
    // Contributions/energy are accumulated by the kernel
    // retreive them and add them to the view
    const double lennardJones =
        computer
            .template getContribution<simd::LJHKernelVecContributions::LJH>();
    const double debyeHueckel =
        computer
            .template getContribution<simd::LJHKernelVecContributions::DH>();

    inView.addContribution(getOffsetContributions(), 0, lennardJones);
    inView.addContribution(getOffsetContributions(), 1, debyeHueckel);
  }

  static std::string GetName() { return "LJH"; }

  std::string type() const override { return GetName(); }

  std::vector<std::string> getExtraContributionLabels() const final {
    std::vector<std::string> labels;
    labels.emplace_back("LJH-lennardJones");
    labels.emplace_back("LJH-debyeHueckel");
    return labels;
  }
};
}  // namespace pairkernels

#endif
