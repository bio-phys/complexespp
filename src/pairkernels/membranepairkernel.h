#ifndef MEMBRANEPAIRKERNEL_H
#define MEMBRANEPAIRKERNEL_H

#include <memory>

#include "abstractpairkernel.h"
#include "constants.h"
#include "domains/abstractdomain.h"
#include "domains/membrane.h"
#include "energy/forcefield.h"

namespace pairkernels {

namespace nc = constants::natural;

// Ulmschneider, M. B., Sansom, M. S. P., & Di Nola, A. (2005). Properties of
// integral membrane protein structures: Derivation of an implicit membrane
// potential. Proteins: Structure, Function and Genetics, 59(2), 252–265.
// https://doi.org/10.1002/prot.20334
inline double f(const double z, const std::array<double, 8>& a) {
  return a[0] + a[1] * std::exp(-a[2] * (z - a[3]) * (z - a[3])) +
         a[4] * std::exp(-a[5] * (z - a[6]) * (z - a[6]));
}

inline double ref(const double z,
                  const std::vector<std::array<double, 8>>& ff) {
  double sum = 0;
  for (const auto& a : ff) {
    sum += f(z, a);
  }
  return -std::log(sum);
}

inline double pot_mean_force(const double z,
                             const std::vector<std::array<double, 8>>& ff,
                             const int idx) {
  return -std::log(f(z, ff[idx])) - ref(z, ff) - ff[idx][7];
}

inline double potential(const double z,
                        const std::vector<std::array<double, 8>>& ff,
                        const int idx, const double z0) {
  if (std::abs(z) < z0) {
    return 10000;
  }
  // catch case where log in potential get's infinite
  const auto& a = ff[idx];
  if (a[0] == 0 && a[1] == 0 && a[4] == 0) {
    return 0;
  }
  return pot_mean_force(z, ff, idx);
}

// McLaughlin, S. (1989). The electrostatic properties of membranes. Annual
// Review of Biophysics and Biophysical Chemistry, 18, 113–136.
// https://doi.org/10.1146/annurev.biophys.18.1.113
inline double guoy_chapman(const double z, const double ePSI0,
                           const double debyeLength) {
  const auto exp_PSI0 = std::exp(ePSI0 / 2 / constants::units::energy);
  const auto alpha = (exp_PSI0 - 1) / (exp_PSI0 + 1);
  const auto k = 1 / debyeLength;
  const auto aa = alpha * std::exp(-k * z);
  return 2 * constants::units::energy / nc::elementaryCharge *
         std::log((1 + aa) / (1 - aa));
}

class MembranePairKernel : public AbstractPairKernel {
  static std::pair<double, double> coreCompute(
      const domains::AbstractDomain& inMembrane,
      const std::pair<int, int>& inInterval, const domains::AbstractDomain& dom,
      const energy::ForceField& forcefield, const util::rvec& box) {
    const auto memFF = forcefield.membrane();
    const auto xyz = dom.xyz();
    const auto beads = dom.beads();
    const auto membrane =
        reinterpret_cast<const domains::Membrane*>(&inMembrane);

    // values in the comment give defaults used in original complexes paper
    const auto z0 = membrane->z0();        // 20.0
    const auto ePSI0 = membrane->ePSI0();  /// -30 [mV]

    double en = 0;
    double en_el = 0;
    for (auto i = inInterval.first; i < inInterval.second; ++i) {
      const auto bead = util::rvec(xyz(i, 0), xyz(i, 1), xyz(i, 2));
      const auto dists = membrane->distance(bead, box);
      const auto charge = dom.charges()[i];
      for (const auto d : dists) {
        en += potential(d, memFF, beads[i], z0);
        en_el += charge * guoy_chapman(d, ePSI0, forcefield.debyeLength());
      }
    }
    return {en, en_el};
  }

 public:
  static const int MembraneNbContributions = 2;

  static std::unique_ptr<AbstractPairKernel> BuildMembranePairKernel() {
    return std::make_unique<MembranePairKernel>();
  }

  explicit MembranePairKernel() : AbstractPairKernel(MembraneNbContributions) {}

  void compute(energy::EnergyMatrixBuffer<>::Accesser inView,
               const domains::AbstractDomain& inDom1,
               const domains::AbstractDomain& inDom2, const util::rvec& box,
               const energy::ForceField& forcefield,
               const std::pair<int, int>& inInterval1,
               const std::pair<int, int>& inInterval2) const final {
    UNUSED(forcefield);

    if (inDom1.isMembrane() && inDom2.isMembrane()) {
      throw std::runtime_error(fmt::format(
          "Error both domains cannot be of type membrane,"
          "ids {} and {}",
          inDom1.id(), inDom2.id()));
    }

    if (!(inDom1.isMembrane() || inDom2.isMembrane())) {
      throw std::runtime_error(fmt::format(
          "At least one domain should be of type membrane,"
          "ids {} and {} (original type {} and {}",
          inDom1.id(), inDom2.id(), inDom1.type(), inDom2.type()));
    }

    if (inDom1.isMembrane()) {
      const auto energy =
          coreCompute(inDom1, inInterval2, inDom2, forcefield, box);
      inView.addContribution(getOffsetContributions(), 0, energy.first);
      inView.addContribution(getOffsetContributions(), 1, energy.second);
    } else {
      const auto energy =
          coreCompute(inDom2, inInterval1, inDom1, forcefield, box);
      inView.addContribution(getOffsetContributions(), 0, energy.first);
      inView.addContribution(getOffsetContributions(), 1, energy.second);
    }
  }

  static std::string GetName() { return "Membrane"; }

  std::string type() const override { return GetName(); }

  std::vector<std::string> getExtraContributionLabels() const final {
    std::vector<std::string> labels;
    labels.emplace_back("Membrane_pot");
    labels.emplace_back("Membrane_el");
    return labels;
  }
};
}  // namespace pairkernels

#endif
