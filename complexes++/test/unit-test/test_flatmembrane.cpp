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

#include "membrane/flatmembrane.h"

TEST(FLATMEMBRANE, distance) {
  const auto z = 50;

  const auto membrane = membrane::Flat<double>(z);

  const auto bead = util::rvec(0, 0, 0);
  const auto box = util::rvec(100, 100, 100);

  const auto dist = membrane.distance(bead, box);
  EXPECT_EQ(dist.size(), 1);
  EXPECT_DOUBLE_EQ(dist[0], 50);
}
