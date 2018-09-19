// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#include <cstdint>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "util/array.h"
#include "util/random.h"

namespace util {

std::uint32_t fullEntropySeed(std::uint32_t seed) noexcept {
  // libc++ from OSX 10.10 doesn't have a public move constructor.
  std::seed_seq seq{seed};
  auto seeds = std::vector<uint32_t>{0};
  seq.generate(seeds.begin(), seeds.end());
  return seeds.at(0);
}

std::string getRNGState(const RNGEngine& rng) {
  std::stringstream iss;
  iss << rng;
  std::string state;
  std::getline(iss, state);
  return state;
}
void setRNGState(RNGEngine& rng, const std::string& state) {
  std::istringstream iss(state);
  iss >> rng;
}
}  // namespace util
