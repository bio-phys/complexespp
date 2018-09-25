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
#ifndef REPULSIVEKERNELVEC_H
#define REPULSIVEKERNELVEC_H

#include <vector>

#include "abstractkernelvec.h"
#include "constants.h"
#include "energy/energy.h"
#include "energy/forcefield.h"
#include "energy/pairparameter.h"
#include "util/linalg.h"
#include "util/moves.h"
#include "util/random.h"
#include "util/util.h"
#include "vectorization/functionsvec.h"

namespace simd {

enum class RepulsiveKernelVecContributions { En };

template <class VecType>
class RepulsiveKernelVec : public AbstractKernelVec<VecType> {
  using RealType = typename VecType::RealType;
  using MaskType = typename VecType::MaskType;

  // Parameters needed to compute the energy
  const energy::PairParameter<RealType> &m_diameter;

  // Accumulate the contribution of the current kernels
  // into these variables:
  VecType Contribution;

public:
  static const int NbContributions = 1;

  using ContributionIdxType = RepulsiveKernelVecContributions;

  // Keep references to the needed information to compute the energy
  inline RepulsiveKernelVec(const RealType inDebyeLength,
                            const RealType inDielectricConstant,
                            const energy::ForceField &inForcefield)
      : m_diameter(inForcefield.diameter()), Contribution(RealType(0)) {
    UNUSED(inDebyeLength);
    UNUSED(inDielectricConstant);
  }

  inline VecType compute(const VecType r2, const VecType inCharge1,
                         const VecType inCharge2, const int inBeadType1,
                         const int inBeadType2[]) override {
    VecType sigma;
    sigma.setFromIndirectArray(m_diameter.data(inBeadType1), inBeadType2);

    const VecType repul = util::repulsive<VecType>(r2, sigma);
    Contribution += repul;
    return repul;
  }

  template <RepulsiveKernelVecContributions inIdxContribution>
  inline VecType getContribution() const {
    static_assert(
        static_cast<int>(inIdxContribution) < NbContributions,
        "inIdxContribution cannot be greater than the number of contributions");
    if (inIdxContribution == RepulsiveKernelVecContributions::En) {
      return Contribution;
    }
  }
};
} // namespace simd

#endif // REPULSIVEKERNELVEC_H
