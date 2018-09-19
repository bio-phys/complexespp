// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
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
  const energy::PairParameter<RealType>& m_diameter;

  // Accumulate the contribution of the current kernels
  // into these variables:
  VecType Contribution;

 public:
  static const int NbContributions = 1;

  using ContributionIdxType = RepulsiveKernelVecContributions;

  // Keep references to the needed information to compute the energy
  inline RepulsiveKernelVec(const RealType inDebyeLength,
                            const RealType inDielectricConstant,
                            const energy::ForceField& inForcefield)
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
}  // namespace simd

#endif  // REPULSIVEKERNELVEC_H
