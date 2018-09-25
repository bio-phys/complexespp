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

#include "mc/accept_func.h"
#include "util/util.h"

namespace mc {
double glauberAccept(const double deltaE, const double beta) noexcept {
  return 1. / (1 + std::exp(beta * deltaE));
}

double metropolisAccept(const double deltaE, const double beta) noexcept {
  return deltaE <= 0 ? 1 : std::exp(-beta * deltaE);
}

double dynamicAccept(const double deltaE, const double beta) noexcept {
  const double us = 1.13124207139993;
  const double u = beta * deltaE;
  if (deltaE < -us) {
    return 1;
  } else if (deltaE > us) {
    return std::exp(-u);
  } else {
    return std::exp(-u * .5) * (3. / 5 - u * u / 40);
  }
}

double alwaysAccept(const double deltaE, const double beta) noexcept {
  UNUSED(deltaE);
  UNUSED(beta);
  return 1;
}

} // namespace mc
