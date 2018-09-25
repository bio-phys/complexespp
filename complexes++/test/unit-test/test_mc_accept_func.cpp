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

#include "mc/accept_func.h"

TEST(MC_UTIL, glauberAccept) {
  EXPECT_NEAR(0.5, mc::glauberAccept(0, 1), 1e-12);
  EXPECT_NEAR(0, mc::glauberAccept(100, 1), 1e-12);
  EXPECT_NEAR(1, mc::glauberAccept(-100, 1), 1e-12);
}

TEST(MC_UTIL, metropolisAccept) {
  EXPECT_NEAR(1, mc::metropolisAccept(0, 1), 1e-12);
  EXPECT_NEAR(0, mc::metropolisAccept(100, 1), 1e-12);
  EXPECT_NEAR(1, mc::metropolisAccept(-100, 1), 1e-12);
}

TEST(MC_UTIL, alwaysAccept) {
  EXPECT_NEAR(1, mc::alwaysAccept(0, 1), 1e-12);
  EXPECT_NEAR(1, mc::alwaysAccept(100, 1), 1e-12);
  EXPECT_NEAR(1, mc::alwaysAccept(-100, 1), 1e-12);
}
