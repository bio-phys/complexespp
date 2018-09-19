// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#include "complexes_test.h"
#include "gtest/gtest.h"

#include "domains/membrane.h"
#include "membrane/flatmembrane.h"

TEST(DOMAINS, membrane_copy) {
  const auto z = 0.0;
  auto mem1 = domains::Membrane(
      "membi", std::make_unique<membrane::Flat<double>>(z), 13, 0, 0, 42,
      domains::Bead(1), domains::BeadChainID("A", 1));
  const auto mem2 = mem1.copy();
  EXPECT_STREQ(mem1.name().c_str(), mem2->name().c_str());
  EXPECT_EQ(mem1.id(), mem2->id());
  EXPECT_EQ(mem1.typeId(), mem2->typeId());
}
