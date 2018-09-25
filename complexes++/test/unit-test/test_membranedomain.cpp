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
