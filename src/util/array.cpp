// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#include <array>
#include <cstdint>
#include <fmt/format.h>
#include <iostream>
#include <vector>

#include "util/array.h"

namespace std {
std::ostream &operator<<(std::ostream &out, const util::rvec &v) {
  out << fmt::format("{:12.8f},{:12.8f},{:12.8f}", v[0], v[1], v[2]);
  return out;
}  // namespace std
}
