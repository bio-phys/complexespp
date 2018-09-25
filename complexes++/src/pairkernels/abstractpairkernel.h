#ifndef ABSTRACTPAIRKERNEL_H
#define ABSTRACTPAIRKERNEL_H

#include <string>
#include <vector>

#include "domains/abstractdomain.h"
#include "energy/energymatrix.h"

namespace pairkernels {

class AbstractPairKernel {
  const int m_nbContributions;
  int m_offsetContributions;

 public:
  explicit AbstractPairKernel(const int inNbContributions,
                              const int inOffsetContributions = 0)
      : m_nbContributions(inNbContributions),
        m_offsetContributions(inOffsetContributions) {}

  virtual ~AbstractPairKernel() {}

  /** This function computes the energy for the beads
   * in given intervals (inInterval1, inInterval2)
   * between the firsst AbstractDomain and the second AbstractDomain.
   * Both domains might be the same, which in this case means
   * that we ask for the inside energy between the beads.
   * This function is maybe called several times.
   */
  virtual void compute(energy::EnergyMatrixBuffer<>::Accesser inView,
                       const domains::AbstractDomain& inDom1,
                       const domains::AbstractDomain& inDom2,
                       const util::rvec& box,
                       const energy::ForceField& forcefield,
                       const std::pair<int, int>& inInterval1,
                       const std::pair<int, int>& inInterval2) const = 0;

  /** This function computes the energy for all the beads
   * between the first AbstractDomain and the second AbstractDomain.
   * Both domains might be the same, which in this case means
   * that we ask for the inside energy between the beads.
   * This function is called only once, and if it is called then
   * energyForIntervals(other) will not be called.
   * The expected result is equal to:
   * @code energyForIntervals(params, {0, inDom1.nbAtoms()}, {0,
   * inDom2.nbAtoms()});
   */
  template <class ViewType>
  void compute(ViewType&& inView, const domains::AbstractDomain& inDom1,
               const domains::AbstractDomain& inDom2, const util::rvec& box,
               const energy::ForceField& forcefield) const {
    return compute(std::forward<ViewType>(inView), inDom1, inDom2, box,
                   forcefield, {0, inDom1.nBeads()}, {0, inDom2.nBeads()});
  }

  int getNbExtraContributions() const { return m_nbContributions; }

  int getOffsetContributions() const { return m_offsetContributions; }

  void setOffsetContributions(const int inOffsetContributions) {
    m_offsetContributions = inOffsetContributions;
  }

  virtual std::vector<std::string> getExtraContributionLabels() const {
    return std::vector<std::string>();
  }

  virtual std::string type() const = 0;
};
}  // namespace pairkernels

#endif
