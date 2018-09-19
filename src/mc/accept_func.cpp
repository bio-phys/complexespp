// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
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

}  // namespace mc
