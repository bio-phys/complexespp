// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#ifndef CODENSEGRIDCONTAINER_H
#define CODENSEGRIDCONTAINER_H

#include <vector>

#include "cutoffgrid/coabstractgridcontainer.h"
#include "cutoffgrid/cocell.h"
#include "cutoffgrid/codomaincelllink.h"
#include "util/array.h"

namespace cutoffgrid {
/**
 * The CoDenseGridContainer class lets store and retrieve
 * cells in a grid. It allocates a dense grid (3dim array)
 * and so all the cells exist at all time.
 * However, the getNbCells method returns the number of cells
 * that have at least one domain (= non empty cell).
 * The access time is guarantee to be constant.
 */
class CoDenseGridContainer : public CoAbstractGridContainer {
 protected:
  //< The size of the grid (number of cells in each dimension)
  util::vec<int> m_gridSize;
  //< The number of cells that have domains
  size_t m_nbExistingCells;
  //< The grid as an array
  std::vector<CoCell> m_grid;

 public:
  explicit CoDenseGridContainer()
      : m_gridSize({0, 0, 0}), m_nbExistingCells(0) {}

  CoDenseGridContainer(const CoDenseGridContainer&) = delete;
  CoDenseGridContainer& operator=(const CoDenseGridContainer&) = delete;
  CoDenseGridContainer(CoDenseGridContainer&&) = default;
  CoDenseGridContainer& operator=(CoDenseGridContainer&&) = default;

  void reset(const util::vec<int>& inGridSize) final {
    m_nbExistingCells = 0;
    m_gridSize = inGridSize;
    const auto allocatedCells = m_gridSize[0] * m_gridSize[1] * m_gridSize[2];
    m_grid.resize(allocatedCells);

    for (int idxX = 0; idxX < m_gridSize[0]; ++idxX) {
      for (int idxY = 0; idxY < m_gridSize[1]; ++idxY) {
        for (int idxZ = 0; idxZ < m_gridSize[2]; ++idxZ) {
          const auto cellIdx = getIndex(util::vec<int>(idxX, idxY, idxZ));
          m_grid[cellIdx].init(cellIdx, util::vec<int>(idxX, idxY, idxZ));
        }
      }
    }
  }

  size_t getNbCells() const final { return m_nbExistingCells; }

  /** This function return the storage index from the box coordinate */
  long long getIndex(const util::vec<int>& inPos) const final {
    const int IndexX = 0;
    const int IndexY = 1;
    const int IndexZ = 2;
    // We use a linear index in Z major
    return (static_cast<long long>(inPos[IndexX] * m_gridSize[IndexY]) +
            inPos[IndexY]) *
               m_gridSize[IndexZ] +
           inPos[IndexZ];
  }

  int addInterval(const long long inCellIdx, const util::vec<int>& coordinate,
                  const int inDomainId, const int inBeginingOfInterval,
                  const int inNbElementsInInterval,
                  const int inPosInCellList) final {
    CoCell& cell = m_grid[inCellIdx];
    if (cell.isEmpty()) {
      m_nbExistingCells += 1;
    }
    const int insertedPos =
        cell.addInterval(inDomainId, inBeginingOfInterval,
                         inNbElementsInInterval, inPosInCellList);
    return insertedPos;
  }

  RemovedInterval removeInterval(
      const CoDomainCellLink& intervalToRemove) final {
    CoCell& cell = m_grid[intervalToRemove.getCellIndex()];
    // The cell manages the deletion of the interval
    const auto resultingUpdate = cell.removeInterval(intervalToRemove);
    if (cell.isEmpty()) {
      m_nbExistingCells -= 1;
    }
    // We return the deletion result unchanged
    return resultingUpdate;
  }

  const CoCell& getCell(const long long inCellIdx) const final {
    return m_grid[inCellIdx];
  }

  const CoCell& getCell(const util::vec<int>& inPos) const final {
    return getCell(getIndex(inPos));
  }

  /////////////////////////////////////////////////////////////////
  /// Iterator over the cells
  /////////////////////////////////////////////////////////////////

  class ConstIterator {
    const CoCell* currentCell;
    const CoCell* const endCell;

    explicit ConstIterator(const CoCell* inCurrentCell, const CoCell* inEndCell)
        : currentCell(inCurrentCell), endCell(inEndCell) {
      while (currentCell != endCell && currentCell->isEmpty()) {
        ++currentCell;
      }
    }

   public:
    const CoCell& operator*() const { return *currentCell; }

    const CoCell* operator->() const { return currentCell; }

    ConstIterator& operator++() {
      DEBUG_ASSERT(currentCell != endCell,
                   "Cannot inc an iterator equal to end()");
      do {
        ++currentCell;
      } while (currentCell != endCell && currentCell->isEmpty());
      return *this;
    }

    ConstIterator operator++(int) {
      DEBUG_ASSERT(currentCell != endCell,
                   "Cannot inc an iterator equal to end()");
      ConstIterator cpIter(currentCell, endCell);
      ++(*this);
      return cpIter;
    }

    bool operator!=(const ConstIterator& other) const {
      return currentCell != other.currentCell;
    }

    bool operator==(const ConstIterator& other) const {
      return currentCell == other.currentCell;
    }

    friend class CoDenseGridContainer;
  };

  ConstIterator cbegin() const {
    return ConstIterator(&m_grid[0], &m_grid[m_grid.size()]);
  }

  ConstIterator cend() const {
    return ConstIterator(&m_grid[m_grid.size()], &m_grid[m_grid.size()]);
  }
};
}  // namespace cutoffgrid
#endif  // CODENSEGRIDCONTAINER_H
