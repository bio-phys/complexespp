// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#include "gtest/gtest.h"

#include "energy/energymatrix.h"

template <int BlockingSize>
void testallCoreLarge() {
  for (int nbRows = 1; nbRows < 1000; nbRows *= 10) {
    for (int nbCols = 1; nbCols < 1000; nbCols *= 10) {
      energy::EnergyMatrix<int, BlockingSize> array(nbRows, nbCols, 0);

      for (int idxRow = 0; idxRow < nbRows; ++idxRow) {
        for (int idxCol = 0; idxCol < nbCols; ++idxCol) {
          EXPECT_EQ(array.getEnergy(idxRow, idxCol), 0);
        }
      }
    }
  }
}

template <int BlockingSize>
void testallCoreSmall() {
  for (int nbRows = 0; nbRows < 64; nbRows += 8) {
    for (int nbCols = 0; nbCols < 64; nbCols += 32) {
      energy::EnergyMatrix<int, BlockingSize> array(nbRows, nbCols, 0);

      for (int idxRow = 0; idxRow < nbRows; ++idxRow) {
        for (int idxCol = 0; idxCol < nbCols; ++idxCol) {
          EXPECT_EQ(array.getEnergy(idxRow, idxCol), 0);
        }
      }
    }
  }
}

template <int BlockingSize>
void testallCoreValues(const int nbValues) {
  for (int nbRows = 0; nbRows < 64; nbRows += 8) {
    for (int nbCols = 0; nbCols < 64; nbCols += 32) {
      energy::EnergyMatrix<int, BlockingSize> array(nbRows, nbCols, nbValues);

      for (int idxRow = 0; idxRow < nbRows; ++idxRow) {
        for (int idxCol = 0; idxCol < nbCols; ++idxCol) {
          EXPECT_FALSE(array.getEnergy(idxRow, idxCol));
          for (int idxValue = 0; idxValue < nbValues; ++idxValue) {
            EXPECT_FALSE(array.getContribution(idxRow, idxCol, idxValue));
          }
        }
      }

      for (int idxRow = 0; idxRow < nbRows; ++idxRow) {
        for (int idxCol = 0; idxCol < nbCols; ++idxCol) {
          for (int idxValue = 0; idxValue < nbValues; ++idxValue) {
            array.addContribution(idxRow, idxCol, idxValue, 1);

            for (int idxRow2 = 0; idxRow2 < nbRows; ++idxRow2) {
              for (int idxCol2 = 0; idxCol2 < nbCols; ++idxCol2) {
                for (int idxValue2 = 0; idxValue2 < nbValues; ++idxValue2) {
                  if (idxRow2 != idxRow || idxCol2 != idxCol ||
                      idxValue2 != idxValue) {
                    EXPECT_FALSE(
                        array.getContribution(idxRow2, idxCol2, idxValue2));
                  } else {
                    EXPECT_TRUE(
                        array.getContribution(idxRow2, idxCol2, idxValue2));
                  }
                }
                if (idxRow2 != idxRow || idxCol2 != idxCol) {
                  EXPECT_FALSE(array.getEnergy(idxRow2, idxCol2));
                } else {
                  EXPECT_TRUE(array.getEnergy(idxRow2, idxCol2));
                }
                EXPECT_FALSE(array.getEnergyConnections(idxRow2, idxCol2));
              }
            }

            array.addContribution(idxRow, idxCol, idxValue, -1);
          }
        }
      }
    }
  }
}

template <int BlockingSize>
void testallCoreReplace(const int nbValues) {
  for (int nbRows = 0; nbRows < 64; nbRows += 8) {
    const int nbCols = nbRows;
    energy::EnergyMatrix<int, BlockingSize> array(nbRows, nbCols, nbValues);

    for (int idxRow = 0; idxRow < nbRows; ++idxRow) {
      for (int idxCol = 0; idxCol < nbCols; ++idxCol) {
        EXPECT_FALSE(array.getEnergy(idxRow, idxCol));
        for (int idxValue = 0; idxValue < nbValues; ++idxValue) {
          EXPECT_FALSE(array.getContribution(idxRow, idxCol, idxValue));
        }
      }
    }

    {
      energy::EnergyMatrix<int, BlockingSize> newValues(nbRows, 1, nbValues);

      newValues.setAll(1);

      for (int idxRow = 0; idxRow < nbRows; ++idxRow) {
        const int idxCol = 0;
        EXPECT_TRUE(newValues.getEnergy(idxRow, idxCol));
        for (int idxValue = 0; idxValue < nbValues; ++idxValue) {
          EXPECT_TRUE(newValues.getContribution(idxRow, idxCol, idxValue));
        }
      }

      for (int idxRow = 0; idxRow < nbRows; ++idxRow) {
        array.replaceRow(idxRow, newValues);
        for (int idxRow2 = 0; idxRow2 < nbRows; ++idxRow2) {
          for (int idxCol2 = 0; idxCol2 < nbCols; ++idxCol2) {
            if (idxRow2 == idxRow) {
              EXPECT_TRUE(array.getEnergy(idxRow2, idxCol2));
              for (int idxValue2 = 0; idxValue2 < nbValues; ++idxValue2) {
                EXPECT_TRUE(array.getContribution(idxRow2, idxCol2, idxValue2));
              }
            } else {
              EXPECT_FALSE(array.getEnergy(idxRow2, idxCol2));
              for (int idxValue2 = 0; idxValue2 < nbValues; ++idxValue2) {
                EXPECT_FALSE(
                    array.getContribution(idxRow2, idxCol2, idxValue2));
              }
            }
          }
        }
        array.reset();
      }

      for (int idxCol = 0; idxCol < nbCols; ++idxCol) {
        array.replaceCol(idxCol, newValues);
        for (int idxRow2 = 0; idxRow2 < nbRows; ++idxRow2) {
          for (int idxCol2 = 0; idxCol2 < nbCols; ++idxCol2) {
            if (idxCol2 == idxCol) {
              EXPECT_TRUE(array.getEnergy(idxRow2, idxCol2));
              for (int idxValue2 = 0; idxValue2 < nbValues; ++idxValue2) {
                EXPECT_TRUE(array.getContribution(idxRow2, idxCol2, idxValue2));
              }
            } else {
              EXPECT_FALSE(array.getEnergy(idxRow2, idxCol2));
              for (int idxValue2 = 0; idxValue2 < nbValues; ++idxValue2) {
                EXPECT_FALSE(
                    array.getContribution(idxRow2, idxCol2, idxValue2));
              }
            }
          }
        }
        array.reset();
      }
    }
  }
}

TEST(BLOCKEDARRAY_TEST, testallSmall) {
  testallCoreSmall<1>();
  testallCoreSmall<3>();
  testallCoreSmall<16>();
  testallCoreSmall<32>();
  testallCoreSmall<128>();
}

TEST(BLOCKEDARRAY_TEST, testallLarge) {
  testallCoreLarge<1>();
  testallCoreLarge<3>();
  testallCoreLarge<16>();
  testallCoreLarge<32>();
  testallCoreLarge<128>();
}

TEST(BLOCKEDARRAY_TEST, testallValues) {
  for (int idxValues = 1; idxValues < 5; ++idxValues) {
    testallCoreValues<1>(idxValues);
    testallCoreValues<3>(idxValues);
    testallCoreValues<16>(idxValues);
    testallCoreValues<32>(idxValues);
    testallCoreValues<128>(idxValues);
  }
}

TEST(BLOCKEDARRAY_TEST, testReplace) {
  for (int idxValues = 1; idxValues < 5; ++idxValues) {
    testallCoreReplace<1>(idxValues);
    testallCoreReplace<3>(idxValues);
    testallCoreReplace<16>(idxValues);
    testallCoreReplace<32>(idxValues);
    testallCoreReplace<128>(idxValues);
  }
}
