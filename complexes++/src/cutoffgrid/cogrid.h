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
#ifndef COGRID_HPP
#define COGRID_HPP

#include <unordered_map>
#include <vector>

#include "cutoffgrid/codensegridcontainer.h"

#include "cutoffgrid/cocell.h"
#include "cutoffgrid/codomaincelllink.h"
#include "domains/abstractdomain.h"
#include "util/array.h"
#include "util/log.h"
#include "util/util.h"

namespace cutoffgrid {
/**
 * The CoGrid class represents a m_grid from the algorithm view.
 * It lets add or remove or move domains and manages the lists.
 * It proposes methods to access a list of CoDomainCellLink for each cell
 * and provides the functions to access the neighbors of a cell.
 * The underlying container for the m_grid can be templatized
 * (GridContainerClass)
 * but we use a dense m_grid by default.
 */
template <class RealType, class GridContainerClass> class CoGrid {
protected:
  //< To avoid to use 3
  static const int m_Dimension = 3;
  //< To avoid to use 0
  static const int m_IndexX = 0;
  //< To avoid to use 1
  static const int m_IndexY = 1;
  //< To avoid to use 2
  static const int m_IndexZ = 2;
  //< The number of neighbors in 3D
  static const int m_NbNeighbors = 26;

  //< The cutoff radius
  RealType m_coRadius;
  //< The size of simulation box
  util::vec<RealType> m_boxSize;
  //< The rounded radius per dimension
  util::vec<RealType> m_coRadiusByDimension;
  //< The number of cells/boxes in each dimension
  util::vec<int> m_gridSize;
  //< The m_grid container
  GridContainerClass m_grid;
  //< The number of domains
  int m_nbDomains;
  //< The CoDomainCellLink list per domain
  std::vector<std::vector<CoDomainCellLink>> m_domainsCells;

public:
  using ConstIterator = typename GridContainerClass::ConstIterator;

  explicit CoGrid(const RealType inCoRadius,
                  const util::vec<RealType> &inBoxSize, const int inNbDomains)
      : m_coRadius(inCoRadius), m_boxSize(inBoxSize), m_nbDomains(inNbDomains) {
    for (int idxDim = 0; idxDim < m_Dimension; ++idxDim) {
      // We must find the minimum radius r' which is a factor of the simulation
      // box:
      // r <= r' , r' * k = SIZE
      m_gridSize[idxDim] =
          std::max(static_cast<int>(m_boxSize[idxDim] / m_coRadius), 1);
      if (m_gridSize[idxDim] < 3) {
        util::Log(
            "The radius of the cutoff should be less than a 1/3 of the box "
            "width.\n"
            "It seems that from the given configuration the Full method is "
            "more appropriate.\n");
      }
      // Then we get the number of boxes/cells in this dimension
      m_coRadiusByDimension[idxDim] =
          m_boxSize[idxDim] / static_cast<RealType>(m_gridSize[idxDim]);
    }

    m_grid.reset(m_gridSize);
    m_domainsCells = std::vector<std::vector<CoDomainCellLink>>(m_nbDomains);
  }

  CoGrid(const CoGrid &) = delete;
  CoGrid &operator=(const CoGrid &) = delete;

  CoGrid(CoGrid &&) = default;
  CoGrid &operator=(CoGrid &&) = default;

  util::vec<int> getCoordFromPosition(const RealType inX, const RealType inY,
                                      const RealType inZ) const {
    const int xCoord = static_cast<int>(inX / m_coRadiusByDimension[m_IndexX]);
    const int yCoord = static_cast<int>(inY / m_coRadiusByDimension[m_IndexY]);
    const int zCoord = static_cast<int>(inZ / m_coRadiusByDimension[m_IndexZ]);
    return util::vec<int>(
        (xCoord + m_gridSize[m_IndexX]) % m_gridSize[m_IndexX],
        (yCoord + m_gridSize[m_IndexY]) % m_gridSize[m_IndexY],
        (zCoord + m_gridSize[m_IndexZ]) % m_gridSize[m_IndexZ]);
  }

  void addDomain(const std::unique_ptr<domains::AbstractDomain> &inDomain) {
    std::vector<CoDomainCellLink> &currentDomainCells =
        m_domainsCells[inDomain->id()];
    if (!currentDomainCells.empty()) {
      throw std::invalid_argument(
          "A domain with the same id has already been inserted");
    }
    const auto &elementsPositions = inDomain->xyz();

    long long currentCellIdx = -1;
    util::vec<int> currentCoordinate = {0};
    int startingElementInterval = -1;

    // We iterate over all the beads
    for (int idxElement = 0; idxElement < elementsPositions.rows();
         ++idxElement) {
      // We compute the cell/box for the current AA
      const auto coordinate =
          getCoordFromPosition(elementsPositions(idxElement, m_IndexX),
                               elementsPositions(idxElement, m_IndexY),
                               elementsPositions(idxElement, m_IndexZ));
      // The m_grid convert the coordinate into index
      const long long elementCellIndex = m_grid.getIndex(coordinate);
      // If it is equal to the current interval index, we keep track of it
      if (elementCellIndex != currentCellIdx) {
        // Otherwise save the current interval
        if (currentCellIdx != -1) {
          const int domainInsertPosInList = m_grid.addInterval(
              currentCellIdx, currentCoordinate, inDomain->id(),
              startingElementInterval, idxElement - startingElementInterval,
              static_cast<int>(currentDomainCells.size()));
          currentDomainCells.emplace_back(currentCellIdx,
                                          domainInsertPosInList);
        }
        // Restart an interval
        currentCellIdx = elementCellIndex;
        currentCoordinate = coordinate;
        startingElementInterval = idxElement;
      }
    }
    // Save the ongoing interval
    const int domainInsertPosInList =
        m_grid.addInterval(currentCellIdx, currentCoordinate, inDomain->id(),
                           startingElementInterval,
                           elementsPositions.rows() - startingElementInterval,
                           static_cast<int>(currentDomainCells.size()));
    currentDomainCells.emplace_back(currentCellIdx, domainInsertPosInList);
  }

  void removeDomain(const int inDomainId) {
    std::vector<CoDomainCellLink> &currentDomainCells =
        m_domainsCells[inDomainId];

    // For each link
    for (CoDomainCellLink &currentDomainCell : currentDomainCells) {
      // Ask the m_grid container to remove the corresponding interval
      const auto domainCellToUpdate = m_grid.removeInterval(currentDomainCell);
      // If there is something to propagate (=> a swap has modified an interval)
      if (domainCellToUpdate.modifiedDomainID() != -1) {
        // Propagate the information
        updateDomainCell(domainCellToUpdate, currentDomainCell.getCellIndex());
      }
      currentDomainCell.setInsetedPosInList(-1);
    }

    // Remove the list from the map
    currentDomainCells.clear();
  }

  const std::vector<CoDomainCellLink> &
  getDomainCells(const int inDomainId) const {
    return m_domainsCells[inDomainId];
  }

  void updateDomain(const std::unique_ptr<domains::AbstractDomain> &inDomain) {
    // A dedicated algorithm is not guaranted to be better
    removeDomain(inDomain->id());
    addDomain(inDomain);
  }

  /** This function update the CoDomainCellLink list of the domain
   * inDomainIdToUpdate
   * because an interval have been swap in the cell list.
   * The previous position of the interval was inOldInsertedPos and it is now
   * inNewInsertedPos.
   * The value inBeginingOfInterval ensure to find the correct interval
   */
  void updateDomainCell(const RemovedInterval &inRemInter,
                        const long long inCellIndex) {
    std::vector<CoDomainCellLink> &domainToUpdate =
        m_domainsCells[inRemInter.modifiedDomainID()];
    CoDomainCellLink &currentDomainCell =
        domainToUpdate[inRemInter.posInCellList()];
    DEBUG_ASSERT(currentDomainCell.getCellIndex() == inCellIndex,
                 "The link is not correct (cell indexes are different)");
    DEBUG_ASSERT(currentDomainCell.getInsertPosInList() ==
                     inRemInter.oldPosInList(),
                 "The link is not correct (position in list are different)");
    currentDomainCell.setInsetedPosInList(inRemInter.newIntervalSize());
  }

  /** This function returns the total potential number of cells in the m_grid */
  size_t getNbMaximumCells() const {
    size_t nbCells = 1;
    for (int idxDim = 0; idxDim < m_Dimension; ++idxDim) {
      nbCells *= m_gridSize[idxDim];
    }
    return nbCells;
  }

  /** This function fill the array inNeighbors with the neighbors of inCell.
   * There is a maximum of 26 neighbors (maybe less on the boundary
   * or if we deal with a very small m_grid).
   * The empty cells are excluded.
   */
  int getCellNeighbors(const util::vec<int> &inCellPos,
                       std::array<const CoCell *, 26> *inNeighbors) const {
    static_assert(
        26 == m_NbNeighbors,
        "Update of the dimension imply a different number of neighbors here");

    int cellsCounter = 0;
    for (int idxX = -1; idxX <= 1; ++idxX) {
      const int xIndex = inCellPos[0] + idxX;
      if (xIndex < 0 || m_gridSize[m_IndexX] <= xIndex) {
        continue;
      }

      for (int idxY = -1; idxY <= 1; ++idxY) {
        const int yIndex = inCellPos[1] + idxY;
        if (yIndex < 0 || m_gridSize[m_IndexY] <= yIndex) {
          continue;
        }

        for (int idxZ = -1; idxZ <= 1; ++idxZ) {
          const int zIndex = inCellPos[2] + idxZ;
          if (zIndex < 0 || m_gridSize[m_IndexZ] <= zIndex) {
            continue;
          }

          if (idxX || idxY || idxZ) {
            const CoCell &aNeighbor =
                m_grid.getCell(util::vec<int>(xIndex, yIndex, zIndex));
            // If this cell exists add it
            if (!aNeighbor.isEmpty()) {
              (*inNeighbors)[cellsCounter++] = &aNeighbor;
            }
          }
        }
      }
    }
    return cellsCounter;
  }

  int getCellNeighbors(const CoCell &inCell,
                       std::array<const CoCell *, 26> *inNeighbors) const {
    return getCellNeighbors(
        util::vec<int>(inCell.getX(), inCell.getY(), inCell.getZ()),
        inNeighbors);
  }

  /** This function fill the array inNeighbors with the neighbors of inCell
   * and takes into account the periodic boundary condition.
   * There is a maximum of 26 neighbors.
   * The empty cells are excluded.
   */
  int getPeriodicCellNeighbors(
      const util::vec<int> &inCellPos,
      std::array<const CoCell *, 26> *inNeighbors) const {
    static_assert(
        26 == m_NbNeighbors,
        "Update of the dimension imply a different number of neighbors here");

    // Not on the boundary, then class method
    if (!((!inCellPos[0]) || (!inCellPos[1]) || (!inCellPos[2]) ||
          (inCellPos[0] == m_gridSize[m_IndexX] - 1) ||
          (inCellPos[1] == m_gridSize[m_IndexY] - 1) ||
          (inCellPos[2] == m_gridSize[m_IndexZ] - 1))) {
      return getCellNeighbors(inCellPos, inNeighbors);
    }

    int cellsCounter = 0;
    for (int idxX = -1; idxX <= 1; ++idxX) {
      int xIndex = inCellPos[0] + idxX;
      RealType currentShiftX = 0;
      if (xIndex < 0) {
        xIndex += m_gridSize[m_IndexX];
        currentShiftX = -m_boxSize[m_IndexX];
      } else if (m_gridSize[m_IndexX] <= xIndex) {
        xIndex -= m_gridSize[m_IndexX];
        currentShiftX = m_boxSize[m_IndexX];
      }

      for (int idxY = -1; idxY <= 1; ++idxY) {
        int yIndex = inCellPos[1] + idxY;
        RealType currentShiftY = 0;
        if (yIndex < 0) {
          yIndex += m_gridSize[m_IndexY];
          currentShiftY = -m_boxSize[m_IndexY];
        } else if (m_gridSize[m_IndexY] <= yIndex) {
          yIndex -= m_gridSize[m_IndexY];
          currentShiftY = m_boxSize[m_IndexY];
        }

        for (int idxZ = -1; idxZ <= 1; ++idxZ) {
          int zIndex = inCellPos[2] + idxZ;
          RealType currentShiftZ = 0;
          if (zIndex < 0) {
            zIndex += m_gridSize[m_IndexZ];
            currentShiftZ = -m_boxSize[m_IndexZ];
          } else if (m_gridSize[m_IndexZ] <= zIndex) {
            zIndex -= m_gridSize[m_IndexZ];
            currentShiftZ = m_boxSize[m_IndexZ];
          }

          if (idxX || idxY || idxZ) {
            const CoCell &aNeighbor =
                m_grid.getCell(util::vec<int>(xIndex, yIndex, zIndex));
            if (!aNeighbor.isEmpty()) {
              (*inNeighbors)[cellsCounter++] = &aNeighbor;
            }
          }
        }
      }
    }
    return cellsCounter;
  }

  int getPeriodicCellNeighbors(
      const CoCell &inCell, std::array<const CoCell *, 26> *inNeighbors) const {
    return getPeriodicCellNeighbors(
        util::vec<int>(inCell.getX(), inCell.getY(), inCell.getZ()),
        inNeighbors);
  }

  const CoCell &getCell(const long long inCellIndex) const {
    return m_grid.getCell(inCellIndex);
  }

  const CoCell &getCell(const util::vec<int> &inPos) const {
    return m_grid.getCell(inPos);
  }

  size_t getNbExistingCells() const { return m_grid.getNbCells(); }

  ConstIterator cbegin() const { return m_grid.cbegin(); }

  ConstIterator cend() const { return m_grid.cend(); }

  const util::vec<int> &getGridSize() const { return m_gridSize; }

  const util::vec<RealType> &getRadiusPerDim() const {
    return m_coRadiusByDimension;
  }

  // For range-loop access

  ConstIterator begin() const { return m_grid.cbegin(); }

  ConstIterator end() const { return m_grid.cend(); }
};
} // namespace cutoffgrid
#endif // COGRID_HPP
