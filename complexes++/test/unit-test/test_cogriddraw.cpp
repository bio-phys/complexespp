// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#include <vector>

#include "gtest/gtest.h"

#include "complexes_test.h"
#include "cutoffgrid/cogrid.h"
#include "cutoffgrid/cosparsegridcontainer.h"
#include "domains/rigiddomain.h"
#include "util/linalg.h"

std::unique_ptr<domains::Rigid> dummyRigid(const int id, const int nbeads,
                                           const domains::Bead bead,
                                           const double charge) {
  const auto Dtrans = util::rvec(.2, .2, .2);
  const auto phi = 2.0;

  const auto range = util::arange<int>(nbeads);
  auto beadChainIDs =
      std::vector<domains::BeadChainID>(nbeads, domains::BeadChainID("A", 0));
  std::transform(std::begin(range), std::end(range), std::begin(beadChainIDs),
                 [](const int i) { return domains::BeadChainID("A", i); });

  auto dom = std::make_unique<domains::Rigid>(
      "none", 0, id, std::vector<domains::Bead>(nbeads, domains::Bead(bead)),
      std::vector<double>(nbeads, charge), beadChainIDs,
      domains::Connections(0), Dtrans, phi);

  return dom;
}

void testGridGridSize(const double inCoRadius, const util::rvec& inBoxSize,
                      const util::vec<int>& inExpectedGridSize) {
  cutoffgrid::CoGrid<double, cutoffgrid::CoSparseGridContainer> grid(
      inCoRadius, inBoxSize, 1);
  EXPECT_EQ(grid.getGridSize()[0], inExpectedGridSize[0]);
  EXPECT_EQ(grid.getGridSize()[1], inExpectedGridSize[1]);
  EXPECT_EQ(grid.getGridSize()[2], inExpectedGridSize[2]);

  if (inBoxSize[0] >= inCoRadius) {
    EXPECT_GE(grid.getRadiusPerDim()[0], inCoRadius);
  }
  if (inBoxSize[1] >= inCoRadius) {
    EXPECT_GE(grid.getRadiusPerDim()[1], inCoRadius);
  }
  if (inBoxSize[2] >= inCoRadius) {
    EXPECT_GE(grid.getRadiusPerDim()[2], inCoRadius);
  }
}

TEST(CUTOFF_COGRIDDRAW, testGridBasic) {
  const double coRadius = 10;
  testGridGridSize(coRadius, util::rvec(100, 200, 300),
                   util::vec<int>(10, 20, 30));
  testGridGridSize(coRadius, util::rvec(190, 201, 300.1),
                   util::vec<int>(19, 20, 30));
  testGridGridSize(coRadius, util::rvec(30, 30, 30), util::vec<int>(3, 3, 3));
}

TEST(CUTOFF_COGRIDDRAW, testGridParticles) {
  const double coRadius = 10;
  const auto boxSize = util::rvec(100, 200, 300);
  const auto boxCenter =
      util::rvec(boxSize[0] / 2, boxSize[1] / 2, boxSize[2] / 2);
  const int nbDomains = 100;
  cutoffgrid::CoGrid<double, cutoffgrid::CoSparseGridContainer> grid(
      coRadius, boxSize, nbDomains);

  const auto radiusDim = grid.getRadiusPerDim();

  for (int idxX = 0; idxX < grid.getGridSize()[0]; ++idxX) {
    const double posX = static_cast<double>(idxX) * radiusDim[0] +
                        radiusDim[0] / 2 + boxCenter[0] - boxSize[0] / 2;
    for (int idxY = 0; idxY < grid.getGridSize()[1]; ++idxY) {
      const double posY = static_cast<double>(idxY) * radiusDim[1] +
                          radiusDim[1] / 2 + boxCenter[1] - boxSize[1] / 2;
      for (int idxZ = 0; idxZ < grid.getGridSize()[2]; ++idxZ) {
        const double posZ = static_cast<double>(idxZ) * radiusDim[2] +
                            radiusDim[2] / 2 + boxCenter[2] - boxSize[2] / 2;

        const auto partCoord = grid.getCoordFromPosition(posX, posY, posZ);
        EXPECT_EQ(idxX, partCoord[0]);
        EXPECT_EQ(idxY, partCoord[1]);
        EXPECT_EQ(idxZ, partCoord[2]);
      }
    }
  }
}

TEST(CUTOFF_COGRIDDRAW, testGridParticlesHard) {
  const double coRadius = 10;
  const auto boxSize = util::rvec(200, 200, 200);
  const int nbDomains = 100;
  cutoffgrid::CoGrid<double, cutoffgrid::CoSparseGridContainer> grid(
      coRadius, boxSize, nbDomains);

  auto checkCoordFromPos = [&](const std::array<double, 3>& pos,
                               const std::array<int, 3>& coord) {
    const auto partCoord = grid.getCoordFromPosition(pos[0], pos[1], pos[2]);
    EXPECT_EQ(coord[0], partCoord[0]);
    EXPECT_EQ(coord[1], partCoord[1]);
    EXPECT_EQ(coord[2], partCoord[2]);
  };

  checkCoordFromPos({{5, 5, 5}}, {{0, 0, 0}});
  checkCoordFromPos({{5, 5, 15}}, {{0, 0, 1}});
  checkCoordFromPos({{15, 5, 15}}, {{1, 0, 1}});
  checkCoordFromPos({{15, 15, 5}}, {{1, 1, 0}});
}

TEST(CUTOFF_COGRIDDRAW, testGridListParticles) {
  const double coRadius = 10;
  const auto boxSize = util::rvec(100, 200, 300);
  const auto boxCenter =
      util::rvec(boxSize[0] / 2, boxSize[1] / 2, boxSize[2] / 2);
  const int nbDomains = 1;
  cutoffgrid::CoGrid<double, cutoffgrid::CoSparseGridContainer> grid(
      coRadius, boxSize, nbDomains);

  const auto radiusDim = grid.getRadiusPerDim();

  const auto totalBoxes =
      grid.getGridSize()[0] * grid.getGridSize()[1] * grid.getGridSize()[2];
  util::rArray m_xyz(totalBoxes, 3);

  for (int idxX = 0; idxX < grid.getGridSize()[0]; ++idxX) {
    const double posX = static_cast<double>(idxX) * radiusDim[0] +
                        radiusDim[0] / 2 + boxCenter[0] - boxSize[0] / 2;
    for (int idxY = 0; idxY < grid.getGridSize()[1]; ++idxY) {
      const double posY = static_cast<double>(idxY) * radiusDim[1] +
                          radiusDim[1] / 2 + boxCenter[1] - boxSize[1] / 2;
      for (int idxZ = 0; idxZ < grid.getGridSize()[2]; ++idxZ) {
        const double posZ = static_cast<double>(idxZ) * radiusDim[2] +
                            radiusDim[2] / 2 + boxCenter[2] - boxSize[2] / 2;
        const int partId =
            ((idxX * grid.getGridSize()[1]) + idxY) * grid.getGridSize()[2] +
            idxZ;
        m_xyz(partId, 0) = posX;
        m_xyz(partId, 1) = posY;
        m_xyz(partId, 2) = posZ;
      }
    }
  }
  std::unique_ptr<domains::AbstractDomain> aDomain =
      dummyRigid(0, totalBoxes, domains::Bead(1), 0);
  aDomain->setXyz(m_xyz);
  grid.addDomain(aDomain);

  const std::vector<cutoffgrid::CoDomainCellLink>& links =
      grid.getDomainCells(0);
  EXPECT_EQ(links.size(), totalBoxes);
}

template <typename Grid>
void checkGridCoherency(Grid& grid, domains::Domains& domains,
                        const int fromDomainIdx, const int toDomainIdx,
                        const int nbCellsInDomains = -1) {
  for (int idxDom = fromDomainIdx; idxDom < toDomainIdx; ++idxDom) {
    const std::vector<cutoffgrid::CoDomainCellLink>& links =
        grid.getDomainCells(domains[idxDom]->id());
    if (nbCellsInDomains != -1) {
      EXPECT_EQ(links.size(), nbCellsInDomains);
    }

    for (int idxLk = 0; idxLk < static_cast<int>(links.size()); ++idxLk) {
      const cutoffgrid::CoDomainCellLink& lk = links[idxLk];
      const cutoffgrid::CoInterval& interval =
          grid.getCell(lk.getCellIndex()).getInterval(lk.getInsertPosInList());
      EXPECT_EQ(interval.getBeginingOfInterval(), idxLk);
      EXPECT_EQ(interval.getEndOfInterval(), idxLk + 1);
      EXPECT_EQ(interval.getDomainId(), domains[idxDom]->id());
    }
  }
}

template <typename Grid>
void checkDomainsInOrder(Grid& grid, domains::Domains& domains,
                         const int fromDomainIdx, const int toDomainIdx,
                         const int nbCellsInDomains = -1) {
  for (int idxDom = fromDomainIdx; idxDom < toDomainIdx; ++idxDom) {
    const auto& links = grid.getDomainCells(domains[idxDom]->id());
    if (nbCellsInDomains != -1) {
      EXPECT_EQ(links.size(), nbCellsInDomains);
    }
    for (const auto& lk : links) {
      EXPECT_EQ(lk.getInsertPosInList(), idxDom);
    }
  }
}

template <typename Grid>
void checkDomainsGeneral(Grid& grid, domains::Domains& domains,
                         const int fromDomainIdx, const int toDomainIdx) {
  for (int idxDom = fromDomainIdx; idxDom < toDomainIdx; ++idxDom) {
    EXPECT_EQ(domains[idxDom]->id(), idxDom);
    const auto& links = grid.getDomainCells(domains[idxDom]->id());

    for (const auto& lk : links) {
      const auto& cell = grid.getCell(lk.getCellIndex());
      EXPECT_FALSE(cell.isEmpty());
      const auto& interval = cell.getInterval(lk.getInsertPosInList());
      EXPECT_EQ(interval.getDomainId(), domains[idxDom]->id());
    }
  }
}

template <typename Grid>
void checkGridGeneral(Grid& grid) {
  for (const auto& currentCell : grid) {
    for (int idxInterval = 0; idxInterval < currentCell.getNbDomains();
         ++idxInterval) {
      const auto& inter = currentCell.getInterval(idxInterval);
      const auto& links = grid.getDomainCells(inter.getDomainId());

      EXPECT_NE(links.size(), 0);
      bool validLkFound = false;
      for (const auto& lk : links) {
        if (lk.getCellIndex() == currentCell.getIndex() &&
            lk.getInsertPosInList() == idxInterval) {
          EXPECT_FALSE(validLkFound);
          validLkFound = true;
        }
      }
      EXPECT_TRUE(validLkFound);
    }
  }
}

TEST(CUTOFF_COGRIDDRAW, testGridDomains) {
  const double coRadius = 10;
  const auto boxSize = util::rvec(100, 200, 300);
  const auto boxCenter =
      util::rvec(boxSize[0] / 2, boxSize[1] / 2, boxSize[2] / 2);
  const int nbDomains = 10;
  cutoffgrid::CoGrid<double, cutoffgrid::CoSparseGridContainer> grid(
      coRadius, boxSize, nbDomains);

  const auto radiusDim = grid.getRadiusPerDim();

  const int totalBoxes =
      grid.getGridSize()[0] * grid.getGridSize()[1] * grid.getGridSize()[2];
  util::rArray m_xyz(totalBoxes, 3);

  for (int idxX = 0; idxX < grid.getGridSize()[0]; ++idxX) {
    const double posX = static_cast<double>(idxX) * radiusDim[0] +
                        radiusDim[0] / 2 + boxCenter[0] - boxSize[0] / 2;
    for (int idxY = 0; idxY < grid.getGridSize()[1]; ++idxY) {
      const double posY = static_cast<double>(idxY) * radiusDim[1] +
                          radiusDim[1] / 2 + boxCenter[1] - boxSize[1] / 2;
      for (int idxZ = 0; idxZ < grid.getGridSize()[2]; ++idxZ) {
        const double posZ = static_cast<double>(idxZ) * radiusDim[2] +
                            radiusDim[2] / 2 + boxCenter[2] - boxSize[2] / 2;

        const int partId =
            ((idxX * grid.getGridSize()[1]) + idxY) * grid.getGridSize()[2] +
            idxZ;
        m_xyz(partId, 0) = posX;
        m_xyz(partId, 1) = posY;
        m_xyz(partId, 2) = posZ;
      }
    }
  }

  domains::Domains domains;
  domains.reserve(nbDomains);
  for (int idxDom = 0; idxDom < nbDomains; ++idxDom) {
    domains.emplace_back(dummyRigid(idxDom, totalBoxes, domains::Bead(1), 0));
    domains.back()->setXyz(m_xyz);
    grid.addDomain(domains.back());
  }

  checkGridCoherency(grid, domains, 0, nbDomains, totalBoxes);
  checkDomainsInOrder(grid, domains, 0, nbDomains, totalBoxes);
  checkGridGeneral(grid);

  grid.removeDomain(nbDomains - 1);
  checkGridCoherency(grid, domains, 0, nbDomains - 1, totalBoxes);
  checkDomainsInOrder(grid, domains, 0, nbDomains - 1, totalBoxes);
  checkGridGeneral(grid);

  grid.removeDomain(0);
  checkGridCoherency(grid, domains, 1, nbDomains - 1, totalBoxes);
  checkGridGeneral(grid);

  // Test the update
  {  // Make all the particles in the first leaf
    util::rArray pos = domains[1]->xyz();
    for (int idxPart = 1; idxPart < pos.rows(); ++idxPart) {
      pos(idxPart, 0) = pos(0, 0);
      pos(idxPart, 1) = pos(0, 1);
      pos(idxPart, 2) = pos(0, 2);
    }
    domains[1]->setXyz(pos);
  }
  grid.updateDomain(domains[1]);
  checkDomainsGeneral(grid, domains, 1, nbDomains - 1);
  checkGridGeneral(grid);

  {  // Merge particles two by two
    util::rArray pos = domains[2]->xyz();
    for (int idxPart = 1; idxPart < pos.rows(); idxPart += 2) {
      pos(idxPart, 0) = pos(idxPart - 1, 0);
      pos(idxPart, 1) = pos(idxPart - 1, 1);
      pos(idxPart, 2) = pos(idxPart - 1, 2);
    }
    domains[2]->setXyz(pos);
  }
  grid.updateDomain(domains[2]);
  checkDomainsGeneral(grid, domains, 1, nbDomains - 1);
  checkGridGeneral(grid);
}

TEST(CUTOFF_COGRIDDRAW, testGridDomainsRandom) {
  const double coRadius = 10;
  const auto boxSize = util::rvec(1000, 200, 300);
  const auto boxCenter =
      util::rvec(boxSize[0] / 2, boxSize[1] / 2, boxSize[2] / 2);
  const int nbDomains = 100;

  const int maxDomSize = 1000;
  const int nbRuns = 200;

  cutoffgrid::CoGrid<double, cutoffgrid::CoSparseGridContainer> grid(
      coRadius, boxSize, nbDomains);

  const auto radiusDim = grid.getRadiusPerDim();

  const int totalBoxes =
      grid.getGridSize()[0] * grid.getGridSize()[1] * grid.getGridSize()[2];
  util::rArray m_xyz(totalBoxes, 3);

  for (int idxX = 0; idxX < grid.getGridSize()[0]; ++idxX) {
    const double posX = static_cast<double>(idxX) * radiusDim[0] +
                        radiusDim[0] / 2 + boxCenter[0] - boxSize[0] / 2;
    for (int idxY = 0; idxY < grid.getGridSize()[1]; ++idxY) {
      const double posY = static_cast<double>(idxY) * radiusDim[1] +
                          radiusDim[1] / 2 + boxCenter[1] - boxSize[1] / 2;
      for (int idxZ = 0; idxZ < grid.getGridSize()[2]; ++idxZ) {
        const double posZ = static_cast<double>(idxZ) * radiusDim[2] +
                            radiusDim[2] / 2 + boxCenter[2] - boxSize[2] / 2;

        const int partId =
            ((idxX * grid.getGridSize()[1]) + idxY) * grid.getGridSize()[2] +
            idxZ;
        m_xyz(partId, 0) = posX;
        m_xyz(partId, 1) = posY;
        m_xyz(partId, 2) = posZ;
      }
    }
  }

  std::default_random_engine randEngine(RAND_SEED);
  std::uniform_int_distribution<int> randDomSize(1, maxDomSize);
  std::uniform_int_distribution<int> randPosIdx(0, totalBoxes - 1);

  domains::Domains domains;
  domains.reserve(nbDomains);

  for (int idxDom = 0; idxDom < nbDomains; ++idxDom) {
    const int nbAtomsInDom = randDomSize(randEngine);
    util::rArray current_xyz(nbAtomsInDom, 3);
    for (int idxPart = 0; idxPart < nbAtomsInDom; ++idxPart) {
      const int posPartIdx = randPosIdx(randEngine);
      current_xyz(idxPart, 0) = m_xyz(posPartIdx, 0);
      current_xyz(idxPart, 1) = m_xyz(posPartIdx, 1);
      current_xyz(idxPart, 2) = m_xyz(posPartIdx, 2);
    }
    domains.emplace_back(dummyRigid(idxDom, nbAtomsInDom, domains::Bead(1), 0));
    domains.back()->setXyz(current_xyz);
    grid.addDomain(domains.back());

    checkDomainsGeneral(grid, domains, 0, idxDom + 1);
    checkGridGeneral(grid);
  }

  std::uniform_int_distribution<int> randDom(0, nbDomains - 1);

  for (int idxRun = 0; idxRun < nbRuns; ++idxRun) {
    const int idxDom = randDom(randEngine);
    const int nbAtomsInDom = domains[idxDom]->xyz().rows();
    util::rArray current_xyz(nbAtomsInDom, 3);
    for (int idxPart = 0; idxPart < nbAtomsInDom; ++idxPart) {
      const int posPartIdx = randPosIdx(randEngine);
      current_xyz(idxPart, 0) = m_xyz(posPartIdx, 0);
      current_xyz(idxPart, 1) = m_xyz(posPartIdx, 1);
      current_xyz(idxPart, 2) = m_xyz(posPartIdx, 2);
    }
    domains[idxDom]->setXyz(current_xyz);
    grid.updateDomain(domains[idxDom]);

    checkGridGeneral(grid);
    checkDomainsGeneral(grid, domains, 0, nbDomains);
  }
}

TEST(CUTOFF_COGRIDDRAW, testGridDomainsNeighbors) {
  const double coRadius = 10;
  const auto boxSize = util::rvec(100, 200, 300);
  const auto boxCenter =
      util::rvec(boxSize[0] / 2, boxSize[1] / 2, boxSize[2] / 2);

  const int nbDomains = 10;
  cutoffgrid::CoGrid<double, cutoffgrid::CoSparseGridContainer> grid(
      coRadius, boxSize, nbDomains);

  const auto radiusDim = grid.getRadiusPerDim();

  const int totalBoxes =
      grid.getGridSize()[0] * grid.getGridSize()[1] * grid.getGridSize()[2];
  util::rArray m_xyz(totalBoxes, 3);

  for (int idxX = 0; idxX < grid.getGridSize()[0]; ++idxX) {
    const double posX = static_cast<double>(idxX) * radiusDim[0] +
                        radiusDim[0] / 2 + boxCenter[0] - boxSize[0] / 2;
    for (int idxY = 0; idxY < grid.getGridSize()[1]; ++idxY) {
      const double posY = static_cast<double>(idxY) * radiusDim[1] +
                          radiusDim[1] / 2 + boxCenter[1] - boxSize[1] / 2;
      for (int idxZ = 0; idxZ < grid.getGridSize()[2]; ++idxZ) {
        const double posZ = static_cast<double>(idxZ) * radiusDim[2] +
                            radiusDim[2] / 2 + boxCenter[2] - boxSize[2] / 2;

        const int partId =
            ((idxX * grid.getGridSize()[1]) + idxY) * grid.getGridSize()[2] +
            idxZ;
        m_xyz(partId, 0) = posX;
        m_xyz(partId, 1) = posY;
        m_xyz(partId, 2) = posZ;
      }
    }
  }

  domains::Domains domains;
  domains.reserve(nbDomains);
  for (int idxDom = 0; idxDom < nbDomains; ++idxDom) {
    domains.emplace_back(dummyRigid(idxDom, totalBoxes, domains::Bead(1), 0));
    domains.back()->setXyz(m_xyz);
    grid.addDomain(domains.back());
  }

  int counterCells = 0;
  for (const auto& currentCell : grid) {
    std::array<const cutoffgrid::CoCell*, 26> neighbors;
    const int nbNeighbors =
        grid.getPeriodicCellNeighbors(currentCell, &neighbors);
    EXPECT_EQ(nbNeighbors, 26);

    int neighExists[27] = {0};
    for (int idxNeig = 0; idxNeig < nbNeighbors; ++idxNeig) {
      int xdiff = neighbors[idxNeig]->getX() - currentCell.getX();
      if (xdiff < -1) {
        xdiff += grid.getGridSize()[0];
        EXPECT_EQ(xdiff, 1);
      } else if (xdiff > 1) {
        xdiff -= grid.getGridSize()[0];
        EXPECT_EQ(xdiff, -1);
      }
      int ydiff = neighbors[idxNeig]->getY() - currentCell.getY();
      if (ydiff < -1) {
        ydiff += grid.getGridSize()[1];
        EXPECT_EQ(ydiff, 1);
      } else if (ydiff > 1) {
        ydiff -= grid.getGridSize()[1];
        EXPECT_EQ(ydiff, -1);
      }
      int zdiff = neighbors[idxNeig]->getZ() - currentCell.getZ();
      if (zdiff < -1) {
        zdiff += grid.getGridSize()[2];
        EXPECT_EQ(zdiff, 1);
      } else if (zdiff > 1) {
        zdiff -= grid.getGridSize()[2];
        EXPECT_EQ(zdiff, -1);
      }

      const int neighPos = ((xdiff + 1) * 3 + (ydiff + 1)) * 3 + (zdiff + 1);
      EXPECT_LT(neighPos, 27);
      EXPECT_GE(neighPos, 0);
      EXPECT_EQ(neighExists[neighPos], 0);
      neighExists[neighPos] += 1;
    }

    counterCells += 1;
  }

  EXPECT_EQ(counterCells, grid.getNbExistingCells());
}

TEST(CUTOFF_COGRIDDRAW, testGridDomainsNeighborsNonPer) {
  const double coRadius = 10;
  const auto boxSize = util::rvec(100, 200, 300);
  const auto boxCenter =
      util::rvec(boxSize[0] / 2, boxSize[1] / 2, boxSize[2] / 2);
  const int nbDomains = 10;
  cutoffgrid::CoGrid<double, cutoffgrid::CoSparseGridContainer> grid(
      coRadius, boxSize, nbDomains);

  const auto radiusDim = grid.getRadiusPerDim();

  const int totalBoxes =
      grid.getGridSize()[0] * grid.getGridSize()[1] * grid.getGridSize()[2];
  util::rArray m_xyz(totalBoxes, 3);

  for (int idxX = 0; idxX < grid.getGridSize()[0]; ++idxX) {
    const double posX = static_cast<double>(idxX) * radiusDim[0] +
                        radiusDim[0] / 2 + boxCenter[0] - boxSize[0] / 2;
    for (int idxY = 0; idxY < grid.getGridSize()[1]; ++idxY) {
      const double posY = static_cast<double>(idxY) * radiusDim[1] +
                          radiusDim[1] / 2 + boxCenter[1] - boxSize[1] / 2;
      for (int idxZ = 0; idxZ < grid.getGridSize()[2]; ++idxZ) {
        const double posZ = static_cast<double>(idxZ) * radiusDim[2] +
                            radiusDim[2] / 2 + boxCenter[2] - boxSize[2] / 2;

        const int partId =
            ((idxX * grid.getGridSize()[1]) + idxY) * grid.getGridSize()[2] +
            idxZ;
        m_xyz(partId, 0) = posX;
        m_xyz(partId, 1) = posY;
        m_xyz(partId, 2) = posZ;
      }
    }
  }

  domains::Domains domains;
  domains.reserve(nbDomains);
  for (int idxDom = 0; idxDom < nbDomains; ++idxDom) {
    domains.emplace_back(dummyRigid(idxDom, totalBoxes, domains::Bead(1), 0));
    domains.back()->setXyz(m_xyz);
    grid.addDomain(domains.back());
  }

  int counterCells = 0;
  for (const auto& currentCell : grid) {
    std::array<const cutoffgrid::CoCell*, 26> neighbors;
    const int nbNeighbors = grid.getCellNeighbors(currentCell, &neighbors);
    int nbBoundarySide = 0;
    if (currentCell.getX() == 0 ||
        currentCell.getX() == grid.getGridSize()[0] - 1) {
      nbBoundarySide += 1;
    }
    if (currentCell.getY() == 0 ||
        currentCell.getY() == grid.getGridSize()[1] - 1) {
      nbBoundarySide += 1;
    }
    if (currentCell.getZ() == 0 ||
        currentCell.getZ() == grid.getGridSize()[2] - 1) {
      nbBoundarySide += 1;
    }
    switch (nbBoundarySide) {
      case 3:
        EXPECT_EQ(nbNeighbors, 26 - 9 - 6 - 4);
        break;
      case 2:
        EXPECT_EQ(nbNeighbors, 26 - 9 - 6);
        break;
      case 1:
        EXPECT_EQ(nbNeighbors, 26 - 9);
        break;
      default:
        EXPECT_EQ(nbBoundarySide, 0);
        EXPECT_EQ(nbNeighbors, 26);
    }

    int neighExists[27] = {0};
    for (int idxNeig = 0; idxNeig < nbNeighbors; ++idxNeig) {
      const int xdiff = currentCell.getX() - neighbors[idxNeig]->getX();
      const int ydiff = currentCell.getY() - neighbors[idxNeig]->getY();
      const int zdiff = currentCell.getZ() - neighbors[idxNeig]->getZ();
      const int neighPos = ((xdiff + 1) * 3 + (ydiff + 1)) * 3 + (zdiff + 1);
      EXPECT_EQ(neighExists[neighPos], 0);
      neighExists[neighPos] += 1;
    }

    counterCells += 1;
  }

  EXPECT_EQ(counterCells, grid.getNbExistingCells());
}
