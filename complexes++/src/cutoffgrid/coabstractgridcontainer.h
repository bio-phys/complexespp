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
#ifndef COABSTRACTGRIDCONTAINER_H
#define COABSTRACTGRIDCONTAINER_H

#include "cutoffgrid/cocell.h"
#include "cutoffgrid/codomaincelllink.h"
#include "domains/abstractdomain.h"
#include "util/array.h"

namespace cutoffgrid {
/**
 * This class represents the interface that
 * a grid container must propose.
 * Copy operator are forbiden for now to ensure that
 * no hidden copy are made (because a grid can be enormous).
 * The getIndex provides the abstraction of a 3D indexing into a linear one.
 * This allow to use any desired indexing system.
 */
class CoAbstractGridContainer {
public:
  explicit CoAbstractGridContainer() = default;
  virtual ~CoAbstractGridContainer() {}

  CoAbstractGridContainer(const CoAbstractGridContainer &) = delete;
  CoAbstractGridContainer &operator=(const CoAbstractGridContainer &) = delete;
  CoAbstractGridContainer(CoAbstractGridContainer &&) = default;
  CoAbstractGridContainer &operator=(CoAbstractGridContainer &&) = default;

  virtual void reset(const util::vec<int> &inGridSize) = 0;
  virtual size_t getNbCells() const = 0;
  /** The pair coordinate/index must be unique */
  virtual long long getIndex(const util::vec<int> &pos) const = 0;
  virtual int addInterval(const long long inCellIdx,
                          const util::vec<int> &coordinate,
                          const int inDomainId, const int inBeginingOfInterval,
                          const int inNbElementsInInterval,
                          const int inPosInCellList) = 0;
  virtual RemovedInterval
  removeInterval(const CoDomainCellLink &intervalToRemove) = 0;

  virtual const CoCell &getCell(const long long inCellIdx) const = 0;
  virtual const CoCell &getCell(const util::vec<int> &inPos) const = 0;
};
} // namespace cutoffgrid

#endif // COABSTRACTGRIDCONTAINER_H
