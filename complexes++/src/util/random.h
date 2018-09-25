// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#ifndef UTIL_RANDOM_H
#define UTIL_RANDOM_H

#include <random>
#include <type_traits>

#include "util/array.h"
#include "util/linalg.h"
#include "util/quaternions/quat.h"

class string;

namespace util {

// use mersenne twister as standart RNG.
using RNGEngine = std::mt19937_64;

// use the STL to generate 32bits of good entropy from a small or bad initial
// seed. If you seed a RNG engine please pass the seed through this function
// first to assure that we seed the RNG with high entropy.
std::uint32_t fullEntropySeed(std::uint32_t seed) noexcept;

template <typename Real,
          typename = std::enable_if<std::is_floating_point<Real>::value>>
vec<Real> randomVec(const vec<Real>& scalingFactor, RNGEngine& rng) noexcept {
  auto dist = std::uniform_real_distribution<Real>{-1, 1};
  return vec<Real>(dist(rng) * scalingFactor[0], dist(rng) * scalingFactor[1],
                   dist(rng) * scalingFactor[2]);
}

template <typename Real,
          typename = std::enable_if<std::is_integral<Real>::value>>
std::vector<Real> randomOrderIndices(const Real n, RNGEngine& rng) {
  auto v = arange<Real>(n);
  std::shuffle(v.begin(), v.end(), rng);
  return v;
}

template <typename Real,
          typename = std::enable_if<std::is_floating_point<Real>::value>>
quaternions::Quaternion<Real> randomQuat(const vec<Real> scale,
                                         const vec<Real>& D, const Real dt,
                                         RNGEngine& rng) {
  DEBUG_ASSERT(.5 * dot(D, util::vec<Real>(1, 1, 1)) * dt < 1,
               "D or dt is to big for random quaternion generation");
  const auto ijk = randomVec(scale, rng);
  const auto s = std::sqrt(1 - .5 * dot(D, util::vec<Real>(1, 1, 1)) * dt);
  const auto norm =
      std::sqrt(s * s + ijk[0] * ijk[0] + ijk[1] * ijk[1] + ijk[2] * +ijk[2]);
  auto q = quaternions::Quaternion<Real>(s / norm, ijk[0] / norm, ijk[1] / norm,
                                         ijk[2] / norm);
  return q;
}

template <typename Real,
          typename = std::enable_if<std::is_floating_point<Real>::value>>
util::Matrix<Real> randomRotation(const Real maxPhi, RNGEngine& rng) {
  const auto axis = randomVec<Real>(vec<Real>(1, 1, 1), rng);
  const auto phi = std::uniform_real_distribution<Real>{-maxPhi, maxPhi}(rng);
  return quaternions::Quaternion<Real>(axis, phi).toMat();
}

std::string getRNGState(const RNGEngine& rng);
void setRNGState(RNGEngine& rng, const std::string& state);

}  // namespace util

#endif  // UTIL_RANDOM_H
