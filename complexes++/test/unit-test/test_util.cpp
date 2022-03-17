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
#include "util/util.h"
#include "util/string.h"

#include "complexes_test.h"
#include "gtest/gtest.h"

TEST(UTIL_UTIL_TEST, indexOf) {
  const auto ar = std::array<int, 5>{{1, 2, 3, 4, 5}};
  const auto el_3 = util::indexOf(std::begin(ar), std::end(ar), 3);
  EXPECT_EQ(el_3, 2);

  const auto vec = std::vector<int>{1, 4};
  const auto el_1 = util::indexOf(std::begin(vec), std::end(vec), 1);
  EXPECT_EQ(el_1, 0);

  const auto el_4 = util::indexOf(std::begin(vec), std::end(vec), 4);
  EXPECT_EQ(el_4, 1);

  const auto el_none = util::indexOf(std::begin(vec), std::end(vec), 5);
  EXPECT_EQ(el_none, -1);
}
