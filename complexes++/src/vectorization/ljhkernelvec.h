// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#ifndef LJHKERNELVEC_H
#define LJHKERNELVEC_H

#include <vector>

#include "abstractkernelvec.h"
#include "constants.h"
#include "energy/energy.h"
#include "energy/forcefield.h"
#include "energy/pairparameter.h"
#include "functionsvec.h"
#include "util/linalg.h"
#include "util/moves.h"
#include "util/random.h"
#include "util/util.h"

namespace simd {

enum class LJHKernelVecContributions { LJH, DH };

template <class VecType>
class LJHKernelVec : public AbstractKernelVec<VecType> {
  using RealType = typename VecType::RealType;
  using MaskType = typename VecType::MaskType;

  //< Parameters needed to compute the energy
  const VecType m_debyeLength;
  const VecType m_dielectricConstant;
  const VecType m_debyeHueckelC;

  const energy::PairParameter<RealType>& m_interActionEnergy;
  const energy::PairParameter<RealType>& m_diameter;

  // Accumulate the contribution of the current kernels
  // into these variables:
  VecType ljContribution;
  VecType dhContribution;

 public:
  static const int NbContributions = 2;

  using ContributionIdxType = LJHKernelVecContributions;

  //< Keep references to the needed information to compute the energy
  inline LJHKernelVec(const RealType inDebyeLength,
                      const RealType inDielectricConstant,
                      const energy::ForceField& inForcefield)
      : m_debyeLength(inDebyeLength),
        m_dielectricConstant(inDielectricConstant),
        m_debyeHueckelC(constants::natural::elementaryCharge *
                        constants::natural::elementaryCharge /
                        (4 * M_PI * constants::natural::epsilon_0) /
                        constants::units::angstrom / constants::units::energy),
        m_interActionEnergy(inForcefield.interActionEnergy()),
        m_diameter(inForcefield.diameter()),
        ljContribution(RealType(0)),
        dhContribution(RealType(0)) {}

  /** Return the energy using Lennard Jones and Hueckel formulation */
  inline VecType compute(const VecType r2, const VecType inCharge1,
                         const VecType inCharge2, const int inBeadType1,
                         const int inBeadType2[]) override {
    VecType sigma;
    sigma.setFromIndirectArray(m_diameter.data(inBeadType1), inBeadType2);

    VecType epsilon;
    epsilon.setFromIndirectArray(m_interActionEnergy.data(inBeadType1),
                                 inBeadType2);

    VecType resultsLj = util::lennardJones<VecType>(r2, epsilon, sigma);
    ljContribution += resultsLj;

    VecType resultsDh = util::debyeHueckel<VecType>(
        inCharge1, inCharge2, r2.sqrt(), m_debyeLength, m_dielectricConstant);
    dhContribution += resultsDh;

    return resultsLj + resultsDh;
  }

  template <LJHKernelVecContributions inIdxContribution>
  inline VecType getContribution() const {
    static_assert(
        static_cast<int>(inIdxContribution) < NbContributions,
        "inIdxContribution cannot be greater than the number of contributions");
    if (inIdxContribution == LJHKernelVecContributions::LJH) {
      return ljContribution;
    }
    if (inIdxContribution == LJHKernelVecContributions::DH) {
      return dhContribution;
    }
  }
};
}  // namespace simd

#endif  // LJHKERNELVEC_H
