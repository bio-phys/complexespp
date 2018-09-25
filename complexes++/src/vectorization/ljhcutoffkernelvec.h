// Copyright (c) 2018 the complexes++ development team and contributors
// (see the file AUTHORS for the full list of names)
//
// This file is part of complexes++.
//
// complexes++ is free software: you can redistribute it and/or modify
// it under the terms of the Lesser GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// complexes++ is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with complexes++.  If not, see <https://www.gnu.org/licenses/>
#ifndef LJHCUTOFFKERNELVEC_H
#define LJHCUTOFFKERNELVEC_H

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

enum class LJHCutoffKernelVecContributions { LJH, DH };

template <class VecType>
class LJHCutoffKernelVec : public AbstractKernelVec<VecType> {
  using RealType = typename VecType::RealType;
  using MaskType = typename VecType::MaskType;

  //< Parameters needed to compute the energy
  const VecType m_debyeLength;
  const VecType m_dielectricConstant;
  const VecType m_debyeHueckelC;

  const energy::PairParameter<RealType> &m_interActionEnergy;
  const energy::PairParameter<RealType> &m_diameter;

  const VecType m_rmin_pow2;
  const VecType m_rmax_pow2;

  // Accumulate the contribution of the current kernels
  // into these variables:
  VecType ljContribution;
  VecType dhContribution;

public:
  static const int NbContributions = 2;

  using ContributionIdxType = LJHCutoffKernelVecContributions;

  //< Keep references to the needed information to compute the energy
  inline LJHCutoffKernelVec(const RealType inDebyeLength,
                            const RealType inDielectricConstant,
                            const energy::ForceField &inForcefield)
      : m_debyeLength(inDebyeLength),
        m_dielectricConstant(inDielectricConstant),
        m_debyeHueckelC(constants::natural::elementaryCharge *
                        constants::natural::elementaryCharge /
                        (4 * M_PI * constants::natural::epsilon_0) /
                        constants::units::angstrom / constants::units::energy),
        m_interActionEnergy(inForcefield.interActionEnergy()),
        m_diameter(inForcefield.diameter()), m_rmin_pow2(1.4 * 1.4),
        m_rmax_pow2(1.8 * 1.8), ljContribution(RealType(0)),
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
    VecType fac = util::smooth<VecType>(r2, sigma, m_rmin_pow2, m_rmax_pow2);

    VecType resultsLj = util::lennardJones<VecType>(r2, epsilon, sigma);
    resultsLj *= fac;
    ljContribution += resultsLj;

    VecType resultsDh = util::debyeHueckel<VecType>(
        inCharge1, inCharge2, r2.sqrt(), m_debyeLength, m_dielectricConstant);
    dhContribution += resultsDh;

    return resultsLj + resultsDh;
  }

  template <LJHCutoffKernelVecContributions inIdxContribution>
  inline VecType getContribution() const {
    static_assert(
        static_cast<int>(inIdxContribution) < NbContributions,
        "inIdxContribution cannot be greater than the number of contributions");
    if (inIdxContribution == LJHCutoffKernelVecContributions::LJH) {
      return ljContribution;
    }
    if (inIdxContribution == LJHCutoffKernelVecContributions::DH) {
      return dhContribution;
    }
  }
};
} // namespace simd

#endif
