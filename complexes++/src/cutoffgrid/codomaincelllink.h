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

  CoDomainCellLink(const CoDomainCellLink &) = default;
  CoDomainCellLink &operator=(const CoDomainCellLink &) = default;
  CoDomainCellLink(CoDomainCellLink &) = default;
  CoDomainCellLink &operator=(CoDomainCellLink &&) = default;

  long long getCellIndex() const { return m_cellIndex; }

  void setCellIndex(const long long inCellIndex) { m_cellIndex = inCellIndex; }

  int getInsertPosInList() const { return m_insertedPosInList; }

  void setInsetedPosInList(const int inInsertedPosInList) {
    m_insertedPosInList = inInsertedPosInList;
  }
};
} // namespace cutoffgrid
#endif // DOMAINCELLLINK_H
