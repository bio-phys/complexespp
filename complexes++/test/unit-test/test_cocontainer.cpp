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
#include <unordered_set>

#include "gtest/gtest.h"

#include "cutoffgrid/codensegridcontainer.h"
#include "cutoffgrid/cosparsegridcontainer.h"

template <class ContainerClass>
void testIndexesIntervals(const util::vec<int> &inGridSize) {
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
void testExistingCells(const util::vec<int> &inGridSize) {
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
void testInsertingPos(const util::vec<int> &inGridSize) {
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
