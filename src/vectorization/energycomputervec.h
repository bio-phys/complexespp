// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#ifndef ENERGYCOMPUTERVEC_H
#define ENERGYCOMPUTERVEC_H

#include "domains/abstractdomain.h"
#include "util/linalg.h"
#include "util/moves.h"
#include "util/pbc.h"
#include "util/random.h"
#include "util/util.h"

#include "InastempConfig.h"
#include "SCALAR/InaVecSCALARDouble.hpp"
#include "SCALAR/InaVecSCALARFloat.hpp"
#include "functionsvec.h"

namespace simd {

/**
 * @brief The EnergyComputer class provides an abstraction
 * over the iteration on the domains to compute the energy.
 * The different methods need a kernel which respects
 * the AbstractKernel interface.
 */
template <template <class KVecType> class KernelClass, class VecType>
class EnergyComputerVec {
  using RealType = typename VecType::RealType;

  KernelClass<VecType> kernelv;
  KernelClass<InaVecSCALAR<RealType>> kernel;

  static_assert(
      KernelClass<VecType>::NbContributions ==
          KernelClass<InaVecSCALAR<RealType>>::NbContributions,
      "Scalar and vector kernels must have the same number of contributions");

  // Forbid allocation in the heap (stack only)
  static void* operator new(std::size_t) = delete;
  static void* operator new[](std::size_t) = delete;

 public:
  EnergyComputerVec(const RealType inDebyeLength,
                    const RealType inDielectricConstant,
                    const energy::ForceField& inForcefield)
      : kernelv(inDebyeLength, inDielectricConstant, inForcefield),
        kernel(inDebyeLength, inDielectricConstant, inForcefield) {}

  /** Compute the energy between two domains for the given itervals */
  inline void disctinctDomainsIntervals(const domains::AbstractDomain& domain1,
                                        const domains::AbstractDomain& domain2,
                                        const ::util::rvec& box,
                                        const std::pair<int, int>& interval1,
                                        const std::pair<int, int>& interval2) {
    const auto& positions1 = domain1.xyz();
    const auto& positions2 = domain2.xyz();
    const auto& charges1 = domain1.charges();
    const auto& charges2 = domain2.charges();
    const auto& beads1 = domain1.beads();
    const auto& beads2 = domain2.beads();

    const int lastVecIdx2 =
        interval2.first +
        ((interval2.second - interval2.first) / VecType::VecLength) *
            VecType::VecLength;

    for (auto idx1 = interval1.first; idx1 < interval1.second; ++idx1) {
      const VecType x1v = positions1(idx1, 0);
      const VecType y1v = positions1(idx1, 1);
      const VecType z1v = positions1(idx1, 2);
      const VecType charge1v = charges1[idx1];

      for (auto idx2 = interval2.first; idx2 < lastVecIdx2;
           idx2 += VecType::VecLength) {
        const VecType dxv = simd::util::mic<VecType>(
            x1v - VecType(&positions2(idx2, 0)), VecType(box[0]));
        const VecType dyv = simd::util::mic<VecType>(
            y1v - VecType(&positions2(idx2, 1)), VecType(box[1]));
        const VecType dzv = simd::util::mic<VecType>(
            z1v - VecType(&positions2(idx2, 2)), VecType(box[2]));

        VecType r2v = dxv * dxv + dyv * dyv + dzv * dzv;

        kernelv.compute(r2v, charge1v, VecType(&charges2[idx2]),
                        static_cast<int>(beads1[idx1]),
                        reinterpret_cast<const int*>(&beads2[idx2]));
      }

      for (auto idx2 = lastVecIdx2; idx2 < interval2.second; ++idx2) {
        ::util::rvec d;
        d[0] = positions1(idx1, 0) - positions2(idx2, 0);
        d[1] = positions1(idx1, 1) - positions2(idx2, 1);
        d[2] = positions1(idx1, 2) - positions2(idx2, 2);

        auto r2 = ::util::pbc::DistSquare(d, box);

        static_cast<RealType>(kernel.compute(
            r2, charges1[idx1], charges2[idx2], static_cast<int>(beads1[idx1]),
            reinterpret_cast<const int*>(&beads2[idx2])));
      }
    }
  }

  /** Compute the energy inside a domain for the given itervals */
  inline void selfDomainIntervals(const domains::AbstractDomain& domain,
                                  const std::pair<int, int>& interval1,
                                  const std::pair<int, int>& interval2) {
    const auto& positions = domain.xyz();
    const auto& beads = domain.beads();
    const auto& charges = domain.charges();

    for (auto idx1 = interval1.first; idx1 < interval1.second; ++idx1) {
      const VecType x1v = positions(idx1, 0);
      const VecType y1v = positions(idx1, 1);
      const VecType z1v = positions(idx1, 2);
      const VecType charge1v = charges[idx1];

      const int idx2Starting = std::max(interval2.first, idx1 + 3);
      const int idx2Ending = std::max(idx2Starting, interval2.second);

      const int lastVecIdx2 =
          idx2Starting + ((idx2Ending - idx2Starting) / VecType::VecLength) *
                             VecType::VecLength;

      for (auto idx2 = idx2Starting; idx2 < lastVecIdx2;
           idx2 += VecType::VecLength) {
        const VecType dxv = x1v - VecType(&positions(idx2, 0));
        const VecType dyv = y1v - VecType(&positions(idx2, 1));
        const VecType dzv = z1v - VecType(&positions(idx2, 2));

        const VecType r2v = dxv * dxv + dyv * dyv + dzv * dzv;

        kernelv.compute(r2v, charge1v, VecType(&charges[idx2]),
                        static_cast<int>(beads[idx1]),
                        reinterpret_cast<const int*>(&beads[idx2]));
      }

      for (auto idx2 = lastVecIdx2; idx2 < idx2Ending; ++idx2) {
        const RealType dx = positions(idx1, 0) - positions(idx2, 0);
        const RealType dy = positions(idx1, 1) - positions(idx2, 1);
        const RealType dz = positions(idx1, 2) - positions(idx2, 2);

        const RealType r2 = dx * dx + dy * dy + dz * dz;

        static_cast<RealType>(kernel.compute(
            r2, charges[idx1], charges[idx2], static_cast<int>(beads[idx1]),
            reinterpret_cast<const int*>(&beads[idx2])));
      }
    }
  }

  template <
      typename KernelClass<VecType>::ContributionIdxType inIdxContribution>
  inline RealType getContribution() const {
    // The contributions/energy have been accumulate by the kernels
    // Here we sum the values of the vectorized and scalar kernels
    // and return them as the complete contribution
    const RealType contributionScalar =
        kernel.template getContribution<inIdxContribution>().horizontalSum();
    const RealType contributionVec =
        kernelv.template getContribution<inIdxContribution>().horizontalSum();
    return contributionScalar + contributionVec;
  }
};

template <template <class KVecType> class KernelClass, class RealType>
using EnergyComputerVecBest =
    class EnergyComputerVec<KernelClass, InaVecBestType<RealType>>;
}  // namespace simd

#endif
