// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#include <cmath>

#include "util/array.h"
#include "util/linalg.h"

namespace util {

void distance(util::rvec& d, const int i, const int j,
              const util::rArray& xyz) {
  d[0] = xyz(i, 0) - xyz(j, 0);
  d[1] = xyz(i, 1) - xyz(j, 1);
  d[2] = xyz(i, 2) - xyz(j, 2);
}

double distance(const util::rvec& v1, const util::rvec& v2) {
  return std::sqrt((v2[0] - v1[0]) * (v2[0] - v1[0]) +
                   (v2[1] - v1[1]) * (v2[1] - v1[1]) +
                   (v2[2] - v1[2]) * (v2[2] - v1[2]));
}

void normalize(util::rvec& v) {
  double c = 1 / std::sqrt(dot(v, v));
  v[0] = v[0] * c;
  v[1] = v[1] * c;
  v[2] = v[2] * c;
}

rvec normalizeVector(const rvec& xyz) {
  auto n = rvec(xyz);
  normalize(n);
  return n;
}

}  // namespace util
