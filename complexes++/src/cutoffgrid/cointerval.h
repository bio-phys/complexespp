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
#ifndef COINTERVAL_HPP
#define COINTERVAL_HPP

namespace cutoffgrid {
/**
 * The CoInterval class represent an interval of elements
 * that are included in a cell.
 * For example, if Domain A has elements (abcdefg),
 * and if (bcd) are in a cell c, therefore we create
 * and interval (A, 1, 3), where 1 is the position of b
 * and 3 the number of elements.
 */
class CoInterval {
protected:
  //< The domain related to the current interval
  int m_domainId;
  //< The position of the first element of the element-list
  int m_beginingOfInterval;
  //< The number of elements in the current interval
  int m_nbElementsInInterval;
  //< The position of the corresponding interval in the cell list
  int m_posInCellList;

public:
  explicit CoInterval(const int inDomainId, const int inBeginingOfInterval,
                      const int inNbElementsInInterval,
                      const int inPosInCellList)
      : m_domainId(inDomainId), m_beginingOfInterval(inBeginingOfInterval),
        m_nbElementsInInterval(inNbElementsInInterval),
        m_posInCellList(inPosInCellList) {}

  CoInterval(const CoInterval &) = default;
  CoInterval &operator=(const CoInterval &) = default;

  CoInterval(CoInterval &&) = default;
  CoInterval &operator=(CoInterval &&) = default;

  int getDomainId() const { return m_domainId; }

  int getBeginingOfInterval() const { return m_beginingOfInterval; }

  void setBeginingOfInterval(const int inBeginingOfInterval) {
    m_beginingOfInterval = inBeginingOfInterval;
  }

  int getNbElementsInInterval() const { return m_nbElementsInInterval; }

  void setNbElementsInInterval(const int inNbElementsInInterval) {
    m_nbElementsInInterval = inNbElementsInInterval;
  }

  int getEndOfInterval() const {
    return m_beginingOfInterval + m_nbElementsInInterval;
  }

  void setPosInCellList(const int inPosInCellList) {
    m_posInCellList = inPosInCellList;
  }

  int getPosInCellList() const { return m_posInCellList; }

  std::pair<int, int> getInterval() const {
    return {m_beginingOfInterval,
            m_beginingOfInterval + m_nbElementsInInterval};
  }
};
} // namespace cutoffgrid
#endif // COINTERVAL_HPP
