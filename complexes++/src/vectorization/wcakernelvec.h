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
#ifndef WCAKERNELVEC_H
#define WCAKERNELVEC_H

#include "energy/pairparameter.h"
#include "util/util.h"
#include "vectorization/abstractkernelvec.h"

namespace simd {

enum class WCAKernelVecContributions { WCA, DH };

template <class VecType>
class WCAKernelVec : public AbstractKernelVec<VecType> {
  using RealType = typename VecType::RealType;
  using MaskType = typename VecType::MaskType;

  const energy::PairParameter<RealType> &m_interActionEnergy;
  const energy::PairParameter<RealType> &m_diameter;

  // Accumulate the contribution of the current kernels
  // into these variables:
  VecType ljContribution;

public:
  static const int NbContributions = 1;

  using ContributionIdxType = WCAKernelVecContributions;

  //< Keep references to the needed information to compute the energy
  inline WCAKernelVec(const RealType inDebyeLength,
                      const RealType inDielectricConstant,
                      const energy::ForceField &inForcefield)
      : m_interActionEnergy(inForcefield.interActionEnergy()),
        m_diameter(inForcefield.diameter()), ljContribution(RealType(0)) {
    UNUSED(inDebyeLength);
    UNUSED(inDielectricConstant);
  }

  /** Return the energy using Lennard Jones and Hueckel formulation */
  inline VecType compute(const VecType r2, const VecType inCharge1,
                         const VecType inCharge2, const int inBeadType1,
                         const int inBeadType2[]) override {
    VecType sigma;
    sigma.setFromIndirectArray(m_diameter.data(inBeadType1), inBeadType2);

    VecType epsilon;
    epsilon.setFromIndirectArray(m_interActionEnergy.data(inBeadType1),
                                 inBeadType2);

    VecType resultsLj = util::wca<VecType>(r2, epsilon, sigma);
    ljContribution += resultsLj;

    return resultsLj;
  }

  template <WCAKernelVecContributions inIdxContribution>
  inline VecType getContribution() const {
    static_assert(
        static_cast<int>(inIdxContribution) < NbContributions,
        "inIdxContribution cannot be greater than the number of contributions");
    if (inIdxContribution == WCAKernelVecContributions::WCA) {
      return ljContribution;
    }
  }
};
} // namespace simd

#endif // WCAKERNELVEC_H
