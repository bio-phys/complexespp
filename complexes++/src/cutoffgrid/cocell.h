// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#ifndef COCELL_H
#define COCELL_H

#include <vector>

#include "cutoffgrid/codomaincelllink.h"
#include "cutoffgrid/cointerval.h"
#include "domains/abstractdomain.h"
#include "util/array.h"

namespace cutoffgrid {

class RemovedInterval {
 private:
  int m_modifiedDomainId;
  int m_posInCellList;
  int m_oldPosInList;
  int m_newIntervalSize;

 public:
  explicit RemovedInterval(const int modifiedDomainID_,
                           const int posInCellList_, const int newIntervalSize_,
                           const int oldPosInList_)
      : m_modifiedDomainId(modifiedDomainID_),
        m_posInCellList(posInCellList_),
        m_oldPosInList(oldPosInList_),
        m_newIntervalSize(newIntervalSize_) {}

  int modifiedDomainID() const { return m_modifiedDomainId; }
  int posInCellList() const { return m_posInCellList; }
  int oldPosInList() const { return m_oldPosInList; }
  int newIntervalSize() const { return m_newIntervalSize; }
};

/**
 * The CoCell class contains the intervals for the different
 * domains that pass through a given box.
 * It stores its own index and box coordinate
 * in addition to a list(vector) of intervals.
 * Removing in interval (removeInterval) is done by
 * swapping with the last elements which may result in
 * the need to update the cell-domain corresponding link
 * (the one that is swapped from back to the new empty place).
 */
class CoCell {
 protected:
  //< The global index of the current cell
  long long m_cellIndex;
  //< The box/cell coordinate
  util::vec<int> m_coordinate;
  //< The list of intervals inside the current cell
  std::vector<CoInterval> m_intervals;

 public:
  explicit CoCell()
      : m_cellIndex(-1), m_coordinate({-1, -1, -1}), m_intervals() {}

  void init(const long long inCellIndex, const util::vec<int>& inCoordinate) {
    if (isInitialized()) {
      throw std::runtime_error(
          "You cannot init a cell that has already been initialized");
    }
    m_cellIndex = inCellIndex;
    m_coordinate = inCoordinate;
  }

  CoCell(const CoCell&) = default;
  CoCell& operator=(const CoCell&) = default;

  CoCell(CoCell&&) = default;
  CoCell& operator=(CoCell&&) = default;

  void clear() {
    m_cellIndex = -1;
    m_intervals.clear();
  }

  int getNbDomains() const { return static_cast<int>(m_intervals.size()); }

  const CoInterval& getInterval(const int inInterIdx) const {
    return m_intervals[inInterIdx];
  }

  CoInterval& getInterval(const int inInterIdx) {
    return m_intervals[inInterIdx];
  }

  /** This function just push a new CoInterval in the internal list */
  int addInterval(const int inDomainId, const int inBeginingOfInterval,
                  const int inNbElementsInInterval, const int inPosInCellList) {
    if (!isInitialized()) {
      throw std::runtime_error(
          "You cannot use cell that has not been initialized");
    }
    m_intervals.emplace_back(inDomainId, inBeginingOfInterval,
                             inNbElementsInInterval, inPosInCellList);
    return static_cast<int>(m_intervals.size() - 1);
  }

  /** The remove function perform a swap between the interval to remove and the
   * last one
   * therefore its complexity is constant, but the swap must be propagated.
   * The return is the swapped Domain (or null if the domain to remove was the
   * last one)
   * and the beginning of the interval, the new position and previous position.
   */
  RemovedInterval removeInterval(const CoDomainCellLink& intervalToRemove) {
    // FIXME: have a DEBUG_ASSERT
    if (intervalToRemove.getInsertPosInList() >=
        static_cast<int>(m_intervals.size())) {
      throw std::invalid_argument(
          "The given interval has a wrong inserted position");
    }

    if (intervalToRemove.getInsertPosInList() ==
        static_cast<int>(m_intervals.size() - 1)) {
      m_intervals.pop_back();
      return RemovedInterval(-1, 0, 0, 0);
    } else {
      m_intervals[intervalToRemove.getInsertPosInList()] =
          std::move(m_intervals[m_intervals.size() - 1]);
      m_intervals.pop_back();
      const CoInterval& modifedInterval =
          m_intervals[intervalToRemove.getInsertPosInList()];
      return RemovedInterval(modifedInterval.getDomainId(),
                             modifedInterval.getPosInCellList(),
                             intervalToRemove.getInsertPosInList(),
                             static_cast<int>(m_intervals.size()));
    }
  }

  bool isInitialized() const { return m_cellIndex != -1; }

  bool isEmpty() const { return m_intervals.empty(); }

  long long getIndex() const { return m_cellIndex; }

  int getX() const { return m_coordinate[0]; }

  int getY() const { return m_coordinate[1]; }

  int getZ() const { return m_coordinate[2]; }
};
}  // namespace cutoffgrid
#endif  // COCELL_H
