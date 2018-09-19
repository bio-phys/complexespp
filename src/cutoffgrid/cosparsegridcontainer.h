// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#ifndef COSPARSEGRIDCONTAINER_HPP
#define COSPARSEGRIDCONTAINER_HPP

#include <unordered_map>

#include "cutoffgrid/coabstractgridcontainer.h"
#include "cutoffgrid/cocell.h"
#include "util/array.h"

namespace cutoffgrid {
/**
 * The CoSparseGridContainer class lets store and retrieve
 * cells in a m_grid. It allocates a hashmap
 * such that only the cells that contains a domain
 * are allocated.
 * The access time is nearly constant but the locality is
 * not predictable.
 * The complexity is tied to the default number of buckets
 * and the average number of cells per bucket.
 */
class CoSparseGridContainer : public CoAbstractGridContainer {
 protected:
  //< The size of the grid (the number of cells in each dimension)
  util::vec<int> m_gridSize;
  //< The hashmap to store the cells
  std::unordered_map<long long, CoCell> m_grid;

  CoCell m_empty;

 public:
  explicit CoSparseGridContainer()
      : m_gridSize({0, 0, 0}), m_grid(), m_empty() {}

  CoSparseGridContainer(const CoSparseGridContainer&) = delete;
  CoSparseGridContainer& operator=(const CoSparseGridContainer&) = delete;
  CoSparseGridContainer(CoSparseGridContainer&&) = default;
  CoSparseGridContainer& operator=(CoSparseGridContainer&&) = default;

  void reset(const util::vec<int>& inGridSize) final {
    m_gridSize = inGridSize;
    m_grid.clear();
  }

  size_t getNbCells() const final { return m_grid.size(); }

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
    if (!cell.isInitialized()) {
      cell.init(inCellIdx, coordinate);
    }
    // The cell manages the deletion of the interval
    const int insertedPos =
        cell.addInterval(inDomainId, inBeginingOfInterval,
                         inNbElementsInInterval, inPosInCellList);
    return insertedPos;
  }

  RemovedInterval removeInterval(
      const CoDomainCellLink& intervalToRemove) final {
    const auto cellIdx = intervalToRemove.getCellIndex();
    const auto cellIter = m_grid.find(cellIdx);
    if (cellIter == m_grid.end()) {
      throw std::invalid_argument("The cell cannot be found");
    }

    CoCell& cell = (*cellIter).second;
    const auto resultingUpdate = cell.removeInterval(intervalToRemove);
    if (cell.isEmpty()) {
      // If the deletion makes the cell empty, we remove it from the map
      m_grid.erase(cellIdx);
    }
    return resultingUpdate;
  }

  const CoCell& getCell(const long long inCellIdx) const final {
    const auto cellIter = m_grid.find(inCellIdx);
    if (cellIter != m_grid.cend()) {
      return (*cellIter).second;
    } else {
      return m_empty;
    }
  }

  const CoCell& getCell(const util::vec<int>& inPos) const final {
    return getCell(getIndex(inPos));
  }

  /////////////////////////////////////////////////////////////////
  /// Iterator over the cells
  /////////////////////////////////////////////////////////////////

  class ConstIterator
      : public std::unordered_map<long long, CoCell>::const_iterator {
    using Parent = std::unordered_map<long long, CoCell>::const_iterator;

   public:
    explicit ConstIterator(
        const std::unordered_map<long long, CoCell>::const_iterator& inOther)
        : Parent(inOther) {}

    const CoCell& operator*() const { return (*this)->second; }
  };

  ConstIterator cbegin() const { return ConstIterator(m_grid.cbegin()); }

  ConstIterator cend() const { return ConstIterator(m_grid.cend()); }
};
}  // namespace cutoffgrid
#endif  // COSPARSEGRIDCONTAINER_HPP
