// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#ifndef DOMAINCELLLINK_H
#define DOMAINCELLLINK_H

namespace cutoffgrid {
/**
 * The class CoDomainCellLink represents path
 * that a domain elements are creating over the cells.
 * If a domain D of elements (abcdef) is included
 * by three cells c1 = (ab), c2 = (cde), c3 = (f),
 * then this domain will have three CoDomainCellLink,
 * and each of them will contain the information related to the
 * intervals (begining, nb elemnts):
 * (c1, x) (c2, y) (c3, z)
 *
 * In addition, a CoDomainCellLink contains a m_insertedPosInList
 * index which is used to retrieve the corresponding CoInterval
 * in a cell (x, y and z in the example).
 * Therefore, if an CoInterval is moved in the list the
 * related CoDomainCellLink must be updated too.
 */
class CoDomainCellLink {
 private:
  //< The index of the cell (~id)
  long long m_cellIndex;
  //< The position of the related CoInterval
  int m_insertedPosInList;

 public:
  explicit CoDomainCellLink(long long inCellIndex,
                            const int inInsertedPosInList)
      : m_cellIndex(inCellIndex), m_insertedPosInList(inInsertedPosInList) {}

  CoDomainCellLink(const CoDomainCellLink&) = default;
  CoDomainCellLink& operator=(const CoDomainCellLink&) = default;
  CoDomainCellLink(CoDomainCellLink&) = default;
  CoDomainCellLink& operator=(CoDomainCellLink&&) = default;

  long long getCellIndex() const { return m_cellIndex; }

  void setCellIndex(const long long inCellIndex) { m_cellIndex = inCellIndex; }

  int getInsertPosInList() const { return m_insertedPosInList; }

  void setInsetedPosInList(const int inInsertedPosInList) {
    m_insertedPosInList = inInsertedPosInList;
  }
};
}  // namespace cutoffgrid
#endif  // DOMAINCELLLINK_H
