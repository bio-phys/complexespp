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
#ifndef FUNCTIONSVEC_H
#define FUNCTIONSVEC_H

#include "constants.h"
#include "util/util.h"

namespace simd {

namespace util {

template <class VecType>
inline VecType debyeHueckel(const VecType charge1, const VecType charge2,
                            const VecType r, const VecType debyeLength,
                            const VecType dielectricConstant) noexcept {
  const VecType c1c2 = charge1 * charge2;
  if ((c1c2 == 0.).isAllTrue() == false) {
    // We compute here if we are not going to multiply all by 0
    constexpr double c_scalar = constants::natural::elementaryCharge *
                                constants::natural::elementaryCharge /
                                (4 * M_PI * constants::natural::epsilon_0) /
                                constants::units::angstrom /
                                constants::units::energy;
    static const VecType c = c_scalar;

    return c1c2 * c * (-r / debyeLength).exp() / (r * dielectricConstant);
  }
  return VecType::GetZero();
}

template <class VecType>
inline VecType repulsive(const VecType r2, const VecType sigma) noexcept {
  const VecType ss_ir2 = sigma * sigma / r2;
  const VecType ri6 = ss_ir2 * ss_ir2 * ss_ir2;
  return ri6 * ri6;
}

template <class VecType>
inline VecType smooth(const VecType r2, const VecType sigma,
                      const VecType rmin_pow2,
                      const VecType rmax_pow2) noexcept {
  const VecType sigma2 = sigma * sigma;
  const VecType x = r2 / sigma2;
  const VecType one(1.);
  const VecType zero(0.);

  const VecType diff = rmax_pow2 - rmin_pow2;
  const VecType inRadius = (rmax_pow2 - x) * (rmax_pow2 - x) *
                           (rmax_pow2 + 2 * x - 3 * rmin_pow2) /
                           (diff * diff * diff);

  return VecType::If(x < rmin_pow2)
      .Then(one)
      .ElseIf(x > rmax_pow2)
      .Then(zero)
      .Else(inRadius);
  //    Original was:
  //    if (x < rmin_pow2) {
  //      return 1;
  //    } else if (x > rmax_pow2) {
  //      return 0;
  //    } else {
  //      const VecType diff = rmax_pow2 - rmin_pow2;
  //      return (rmax_pow2 - x) * (rmax_pow2 - x) * (rmax_pow2 + 2 * x - 3 *
  //      rmin_pow2) /
  //             (diff * diff * diff);
  //    }
}

template <class VecType>
inline VecType lennardJones(const VecType r2, const VecType epsilon,
                            const VecType sigma) noexcept {
  static const VecType r2ij(1.2599210498948723);

  const VecType ss_ir2 = sigma * sigma / r2;
  const VecType ri6 = ss_ir2 * ss_ir2 * ss_ir2;
  const VecType ri12 = ri6 * ri6;
  const VecType r2ij0 = r2ij * sigma * sigma;

  // Non optimized is given by:
  //  VecType results = VecType::If(epsilon == 0)
  //                        .Then(.01 * ri12)
  //                        .ElseIf(epsilon < 0)
  //                        .Then(-4.0 * epsilon * (ri12 - ri6))
  //                        .ElseIf(r2 < r2ij0)
  //                        .Then(4.0 * epsilon * (ri12 - ri6 + .5))
  //                        .Else(-4 * epsilon * (ri12 - ri6));

  const VecType coef = -4 * epsilon * (ri12 - ri6);
  VecType results = VecType::If(epsilon == 0)
                        .Then(.01 * ri12)
                        .ElseIf(epsilon > 0 & r2 < r2ij0)
                        .Then(2.0 * epsilon - coef)
                        .Else(coef);

  return results;
}

template <class VecType>
inline VecType wca(const VecType r2, const VecType epsilon,
                   const VecType sigma) noexcept {
  static const VecType r2ij(1.2599210498948723);

  const VecType ss_ir2 = sigma * sigma / r2;
  const VecType ri6 = ss_ir2 * ss_ir2 * ss_ir2;
  const VecType ri12 = ri6 * ri6;
  const VecType r2ij0 = r2ij * sigma * sigma;

  const VecType coef = 4 * epsilon * (ri12 - ri6);
  return VecType::If(r2 < r2ij0).Then(-epsilon).Else(coef);
}

template <class VecType>
inline VecType mic(VecType posdiff, const VecType dim) {
  const VecType halfdim = dim * 0.5;
  typename VecType::MaskType msk = (posdiff > halfdim);
  while (msk.isAllFalse() == false) {
    posdiff -= VecType::IfTrue(msk, dim);
    msk = (posdiff > halfdim);
  }
  msk = (posdiff <= -halfdim);
  while (msk.isAllFalse() == false) {
    posdiff += VecType::IfTrue(msk, dim);
    msk = (posdiff <= -halfdim);
  }
  return posdiff;
}
} // namespace util
} // namespace simd

#endif
