// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------

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
