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

#include "mc/exchangelooper.h"

#include "complexes_test.h"
#include "gtest/gtest.h"

TEST(EXCHANGELOOP, OddEvenNeighborLoop) {
  for (int idxSize = 1; idxSize < 30; ++idxSize) {
    for (int idxSweep = 0; idxSweep <= 40; ++idxSweep) {
      std::vector<bool> tested(idxSize, false);

      int totalLoop = 0;

      for (auto iterLoop : mc::OddEvenNeighborLoop(idxSize, idxSweep)) {
        const int idx1 = iterLoop.first;
        const int idx2 = iterLoop.second;
        EXPECT_FALSE(tested[idx1]);
        EXPECT_FALSE(tested[idx2]);

        tested[idx1] = true;
        tested[idx2] = true;

        EXPECT_EQ(std::abs(idx1 - idx2), 1);

        // Those two lines should do as the commented snipet
        EXPECT_EQ((idx1 < idx2 ? idx1 : idx2) & 1, (idxSweep & 1));
        EXPECT_EQ(std::max(idx1, idx2) & 1, !(idxSweep & 1));
        // Gcc 5 complains:
        // error: assuming signed overflow does not occur when assuming that (X
        // + c) < X is always false
        // if (idxSweep & 1) {
        //  EXPECT_EQ(std::min(idx1, idx2) & 1, 1);
        //  EXPECT_EQ(std::max(idx1, idx2) & 1, 0);
        //} else {
        //  EXPECT_EQ(std::max(idx1, idx2) & 1, 1);
        //  EXPECT_EQ(std::min(idx1, idx2) & 1, 0);
        //}

        totalLoop += 1;
      }

      if (idxSweep & 1) {
        EXPECT_EQ(totalLoop, (idxSize - 1) / 2);
      } else {
        EXPECT_EQ(totalLoop, idxSize / 2);
      }

      if (idxSize > 2) {
        EXPECT_TRUE(tested[idxSize - 2]);
      }
      if (idxSize > 1 && (idxSweep & 1) == (idxSize & 1)) {
        EXPECT_TRUE(tested[idxSize - 1]);
      }
    }
  }
}
