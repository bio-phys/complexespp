// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#include "util/util.h"

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
