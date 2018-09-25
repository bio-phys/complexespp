#ifndef WCAPAIRKERNEL_H
#define WCAPAIRKERNEL_H

#include <memory>
#include <vector>

#include "abstractpairkernel.h"
#include "domains/abstractdomain.h"
#include "energy/forcefield.h"

#include "vectorization/energycomputervec.h"
#include "vectorization/wcakernelvec.h"

namespace pairkernels {

class WCAPairKernel : public AbstractPairKernel {
 public:
  static const int WCANbContributions = 1;

  static std::unique_ptr<AbstractPairKernel> BuildWCAPairKernel() {
    return std::make_unique<WCAPairKernel>();
  }

  explicit WCAPairKernel() : AbstractPairKernel(WCANbContributions) {}

  void compute(energy::EnergyMatrixBuffer<>::Accesser inView,
               const domains::AbstractDomain& inDom1,
               const domains::AbstractDomain& inDom2, const util::rvec& box,
               const energy::ForceField& forceField,
               const std::pair<int, int>& inInterval1,
               const std::pair<int, int>& inInterval2) const final {
    const double debyeLength(forceField.debyeLength());
    const double dielectricConstant(forceField.dielectricConstant());
    simd::EnergyComputerVecBest<simd::WCAKernelVec, double> computer(
        debyeLength, dielectricConstant, forceField);
    computer.disctinctDomainsIntervals(inDom1, inDom2, box, inInterval1,
                                       inInterval2);
    // Contributions/energy are accumulated by the kernel
    // retreive them and add them to the view
    const double lennardJones =
        computer
            .template getContribution<simd::WCAKernelVecContributions::WCA>();

    inView.addContribution(getOffsetContributions(), 0, lennardJones);
  }

  static std::string GetName() { return "WCA"; }
  std::string type() const override { return GetName(); }
  std::vector<std::string> getExtraContributionLabels() const final {
    std::vector<std::string> labels;
    labels.emplace_back("WCA");
    return labels;
  }
};
}  // namespace pairkernels

#endif
