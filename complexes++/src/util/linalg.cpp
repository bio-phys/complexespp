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
#include <cmath>

#include "util/array.h"
#include "util/linalg.h"

namespace util {

void distance(util::rvec &d, const int i, const int j,
              const util::rArray &xyz) {
  d[0] = xyz(i, 0) - xyz(j, 0);
  d[1] = xyz(i, 1) - xyz(j, 1);
  d[2] = xyz(i, 2) - xyz(j, 2);
}

double distance(const util::rvec &v1, const util::rvec &v2) {
  return std::sqrt((v2[0] - v1[0]) * (v2[0] - v1[0]) +
                   (v2[1] - v1[1]) * (v2[1] - v1[1]) +
                   (v2[2] - v1[2]) * (v2[2] - v1[2]));
}

void normalize(util::rvec &v) {
  double c = 1 / std::sqrt(dot(v, v));
  v[0] = v[0] * c;
  v[1] = v[1] * c;
  v[2] = v[2] * c;
}

rvec normalizeVector(const rvec &xyz) {
  auto n = rvec(xyz);
  normalize(n);
  return n;
}

} // namespace util
