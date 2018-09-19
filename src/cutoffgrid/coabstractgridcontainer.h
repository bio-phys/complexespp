// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
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

  CoAbstractGridContainer(const CoAbstractGridContainer&) = delete;
  CoAbstractGridContainer& operator=(const CoAbstractGridContainer&) = delete;
  CoAbstractGridContainer(CoAbstractGridContainer&&) = default;
  CoAbstractGridContainer& operator=(CoAbstractGridContainer&&) = default;

  virtual void reset(const util::vec<int>& inGridSize) = 0;
  virtual size_t getNbCells() const = 0;
  /** The pair coordinate/index must be unique */
  virtual long long getIndex(const util::vec<int>& pos) const = 0;
  virtual int addInterval(const long long inCellIdx,
                          const util::vec<int>& coordinate,
                          const int inDomainId, const int inBeginingOfInterval,
                          const int inNbElementsInInterval,
                          const int inPosInCellList) = 0;
  virtual RemovedInterval removeInterval(
      const CoDomainCellLink& intervalToRemove) = 0;

  virtual const CoCell& getCell(const long long inCellIdx) const = 0;
  virtual const CoCell& getCell(const util::vec<int>& inPos) const = 0;
};
}  // namespace cutoffgrid

#endif  // COABSTRACTGRIDCONTAINER_H
