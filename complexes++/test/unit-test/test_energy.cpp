// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#include "complexes_test.h"
#include "gtest/gtest.h"
#include <cmath>

#include "constants.h"
#include "energy/energy.h"
#include "util/array.h"

namespace nc = constants::natural;

// NOTE ABOUT THE UNITS:
// length [Anstroem]
// energy [kT]

TEST(ENERGYTEST, harmonic) {
  const auto rmin = 1.0;
  const auto k = 2.0;
  EXPECT_DOUBLE_EQ(0, energy::harmonic(1, rmin, k));
  EXPECT_DOUBLE_EQ(1, energy::harmonic(0, rmin, k));
  EXPECT_DOUBLE_EQ(1, energy::harmonic(2, rmin, k));
}
