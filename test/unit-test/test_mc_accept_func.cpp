// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
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
