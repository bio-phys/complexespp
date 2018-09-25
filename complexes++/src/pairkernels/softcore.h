
#ifndef SOFTCOREPAIRKERNEL_H
#define SOFTCOREPAIRKERNEL_H

#include <cmath>
#include <memory>
#include <vector>

#include "constants.h"
#include "domains/abstractdomain.h"
#include "energy/forcefield.h"
#include "pairkernels/abstractpairkernel.h"
#include "util/pbc.h"

namespace pairkernels {

inline double softcore(const double r2, const double sigma,
                       const double epsilon, const double alpha) {
  // function minimum
  constexpr auto sr2ij = 1.2599210498948723;
  const auto s = (std::pow(2, 1. / 6) - std::pow(2 - alpha, 1. / 6)) * sigma;
  const auto s6 = std::pow(sigma, 6);
  const auto r = std::sqrt(r2);

  const auto tmp = s6 / (alpha * s6 + std::pow(r - s, 6));
  if (epsilon < 0) {
    return -4.0 * epsilon * (tmp * tmp - tmp);
  } else {
    if (r2 / sigma / sigma < sr2ij) {
      return 4.0 * epsilon * (tmp * tmp - tmp + .5);
    } else {
      return -4.0 * epsilon * (tmp * tmp - tmp);
    }
  }
}

inline double mixingRule(const double l1, const double l2) {
  return l1 * l2 / (l1 + l2);
}

inline double electric(const double q1, const double q2, const double r,
                       const double lambda, const double zeta,
                       const double dielectricConstant) {
  const auto C = constants::natural::elementaryCharge *
                 constants::natural::elementaryCharge /
                 (4 * M_PI * constants::natural::epsilon_0) /
                 constants::units::angstrom / constants::units::energy /
                 dielectricConstant;

  //  const auto C =
  //      1 / (4 * M_PI * dielectricConstant * constants::natural::epsilon_0);
  return q1 * q2 * C * std::erf(r * std::sqrt(lambda)) / r *
         std::exp(-r / zeta);
}

class SoftcorePairKernel : public AbstractPairKernel {
 public:
  static const int SoftcoreNbContributions = 2;

  static std::unique_ptr<AbstractPairKernel> BuildSoftcorePairKernel() {
    return std::make_unique<SoftcorePairKernel>();
  }

  explicit SoftcorePairKernel() : AbstractPairKernel(SoftcoreNbContributions) {}

  void compute(energy::EnergyMatrixBuffer<>::Accesser inView,
               const domains::AbstractDomain& inDom1,
               const domains::AbstractDomain& inDom2, const util::rvec& box,
               const energy::ForceField& forceField,
               const std::pair<int, int>& inInterval1,
               const std::pair<int, int>& inInterval2) const final {
    const auto debyeLength = forceField.debyeLength();
    const auto dielectricConstant = forceField.dielectricConstant();
    const auto alpha = forceField.alpha();

    const auto& beads1 = inDom1.beads();
    const auto& beads2 = inDom2.beads();
    const auto& positions1 = inDom1.xyz();
    const auto& positions2 = inDom2.xyz();
    const auto& charges1 = inDom1.charges();
    const auto& charges2 = inDom2.charges();

    const auto& diameter = forceField.diameter();
    const auto& epsilon = forceField.interActionEnergy();
    const auto& chargeRadii = forceField.chargeRadius();

    auto sc = 0.0;
    auto el = 0.0;
    for (auto idx1 = inInterval1.first; idx1 < inInterval1.second; ++idx1) {
      for (auto idx2 = inInterval2.first; idx2 < inInterval2.second; ++idx2) {
        const auto sigma = diameter(beads1[idx1], beads2[idx2]);
        const auto eps = epsilon(beads1[idx1], beads2[idx2]);

        const auto p1 = util::rvec(positions1(idx1, 0), positions1(idx1, 1),
                                   positions1(idx1, 2));
        const auto p2 = util::rvec(positions2(idx2, 0), positions2(idx2, 1),
                                   positions2(idx2, 2));
        const auto dist =
            util::rvec(p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2]);
        const auto r2 = util::pbc::DistSquare(dist, box);
        sc += softcore(r2, sigma, eps, alpha);

        if (charges1[idx1] != 0) {
          const auto cr =
              mixingRule(chargeRadii[beads1[idx1]], chargeRadii[beads2[idx2]]);
          const auto r = std::sqrt(r2);
          const auto q1 = charges1[idx1];
          const auto q2 = charges2[idx2];
          el += electric(q1, q2, r, cr, debyeLength, dielectricConstant);
        }
      }
    }
    inView.addContribution(getOffsetContributions(), 0, sc);
    inView.addContribution(getOffsetContributions(), 1, el);
  }

  static std::string GetName() { return "Softcore"; }

  std::string type() const override { return GetName(); }

  std::vector<std::string> getExtraContributionLabels() const final {
    std::vector<std::string> labels;
    labels.emplace_back("softcore-pot");
    labels.emplace_back("softcore-el");
    return labels;
  }
};
}  // namespace pairkernels

#endif
