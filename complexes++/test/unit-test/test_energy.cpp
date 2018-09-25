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
