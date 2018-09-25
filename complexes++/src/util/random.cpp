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

std::string getRNGState(const RNGEngine &rng) {
  std::stringstream iss;
  iss << rng;
  std::string state;
  std::getline(iss, state);
  return state;
}
void setRNGState(RNGEngine &rng, const std::string &state) {
  std::istringstream iss(state);
  iss >> rng;
}
} // namespace util
