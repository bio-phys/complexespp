// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
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
