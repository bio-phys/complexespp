// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#include "complexes_test.h"
#include "gtest/gtest.h"

#include "domains/beads.h"

TEST(BEADS, unkown_bead_type) {
  const auto beadTypes = std::vector<std::string>{{"ALA", "VAL", "PRO"}};
  const auto unkown = std::string{"FOO"};
  EXPECT_THROW(domains::findBeadID(unkown, beadTypes), std::runtime_error);
}

TEST(BEADS, find_bead_type) {
  const auto beadTypes = std::vector<std::string>{{"ALA", "VAL", "PRO"}};
  const auto val = std::string{"VAL"};
  EXPECT_EQ(domains::findBeadID(val, beadTypes), 1);
}
