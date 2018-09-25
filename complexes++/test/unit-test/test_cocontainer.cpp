// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#include <unordered_set>

#include "gtest/gtest.h"

#include "cutoffgrid/codensegridcontainer.h"
#include "cutoffgrid/cosparsegridcontainer.h"

template <class ContainerClass>
void testIndexesIntervals(const util::vec<int>& inGridSize) {
  std::unordered_set<long long> testingSet;

  ContainerClass container;
  container.reset(inGridSize);

  for (int idxX = 0; idxX < inGridSize[0]; ++idxX) {
    for (int idxY = 0; idxY < inGridSize[1]; ++idxY) {
      for (int idxZ = 0; idxZ < inGridSize[2]; ++idxZ) {
        const auto index = container.getIndex(util::vec<int>(idxX, idxY, idxZ));
        EXPECT_EQ(testingSet.find(index), testingSet.end());
        testingSet.insert(index);
      }
    }
  }
}

template <class ContainerClass>
void testExistingCells(const util::vec<int>& inGridSize) {
  std::unordered_set<long long> testingSet;

  ContainerClass container;
  container.reset(inGridSize);

  for (int idxX = 0; idxX < inGridSize[0]; ++idxX) {
    for (int idxY = 0; idxY < inGridSize[1]; ++idxY) {
      for (int idxZ = 0; idxZ < inGridSize[2]; ++idxZ) {
        const auto index = container.getIndex(util::vec<int>(idxX, idxY, idxZ));
        EXPECT_TRUE(container.getCell(index).isEmpty());
        container.addInterval(index, util::vec<int>(idxX, idxY, idxZ), 0, 0, 0,
                              0);
        EXPECT_FALSE(container.getCell(index).isEmpty());
      }
    }
  }

  EXPECT_EQ(container.getNbCells(),
            inGridSize[0] * inGridSize[1] * inGridSize[2]);
}

template <class ContainerClass>
void testInsertingPos(const util::vec<int>& inGridSize) {
  std::unordered_set<long long> testingSet;

  ContainerClass container;
  container.reset(inGridSize);

  for (int idxPos = 0; idxPos < 3; ++idxPos) {
    for (int idxX = 0; idxX < inGridSize[0]; ++idxX) {
      for (int idxY = 0; idxY < inGridSize[1]; ++idxY) {
        for (int idxZ = 0; idxZ < inGridSize[2]; ++idxZ) {
          const auto index =
              container.getIndex(util::vec<int>(idxX, idxY, idxZ));
          EXPECT_TRUE(container.getCell(index).isEmpty() || idxPos);
          EXPECT_EQ(container.addInterval(
                        index, util::vec<int>(idxX, idxY, idxZ), 0, 0, 0, 0),
                    idxPos);
          EXPECT_FALSE(container.getCell(index).isEmpty());
        }
      }
    }

    EXPECT_EQ(container.getNbCells(),
              inGridSize[0] * inGridSize[1] * inGridSize[2]);
  }
}

TEST(CUTOFF_TEST, testDenseContainer) {
  testIndexesIntervals<cutoffgrid::CoDenseGridContainer>(
      util::vec<int>(10, 4, 5));
  testIndexesIntervals<cutoffgrid::CoDenseGridContainer>(
      util::vec<int>(1, 1, 5));

  testExistingCells<cutoffgrid::CoDenseGridContainer>(util::vec<int>(10, 4, 5));
  testExistingCells<cutoffgrid::CoDenseGridContainer>(util::vec<int>(1, 1, 5));

  testInsertingPos<cutoffgrid::CoDenseGridContainer>(util::vec<int>(10, 4, 5));
  testInsertingPos<cutoffgrid::CoDenseGridContainer>(util::vec<int>(1, 1, 5));
}

TEST(CUTOFF_TEST, testSparseContainer) {
  testIndexesIntervals<cutoffgrid::CoSparseGridContainer>(
      util::vec<int>(10, 4, 5));
  testIndexesIntervals<cutoffgrid::CoSparseGridContainer>(
      util::vec<int>(1, 1, 5));

  testExistingCells<cutoffgrid::CoSparseGridContainer>(
      util::vec<int>(10, 4, 5));
  testExistingCells<cutoffgrid::CoSparseGridContainer>(util::vec<int>(1, 1, 5));

  testInsertingPos<cutoffgrid::CoSparseGridContainer>(util::vec<int>(10, 4, 5));
  testInsertingPos<cutoffgrid::CoSparseGridContainer>(util::vec<int>(1, 1, 5));
}
