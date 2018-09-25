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
#ifndef ABSTRACTINTERACTIONALGORITHM_HPP
#define ABSTRACTINTERACTIONALGORITHM_HPP

#include "domains/abstractdomain.h"
#include "energy/energymatrix.h"
#include "io/rebuilder.h"
#include "io/serializer.h"
#include "pairkernels/pairkernelmanager.h"

/**
 * The AbstractInteractionAlgorithm algorithm is used to
 * compute the interactions between domains.
 * It fills a matrice in case of all-with-all interactions
 * and a vector in case of one-with-all.
 * When a domain is updated the method updateDomain must be called.
 */
template <class RealType>
class AbstractInteractionAlgorithm
    : public io::RebuilderCore<AbstractInteractionAlgorithm<RealType>,
                               const std::shared_ptr<domains::Domains> &>,
      public io::AbstractSerializable {
protected:
  void serializeCore(io::Serializer &serializer) const {
    serializer.append(type(), "type");
  }

public:
  AbstractInteractionAlgorithm(io::Deserializer &deserializer) {}
  AbstractInteractionAlgorithm() {}
  virtual ~AbstractInteractionAlgorithm() {}

  virtual void updateDomain(const int inDomainId) = 0;

  virtual void resetAllDomains(const util::rvec &newBox) = 0;

  virtual void computeAll(const util::rvec &box,
                          const energy::ForceField &forcefield,
                          const pairkernels::PairKernelManager &inKernels,
                          energy::EnergyMatrix<RealType> &outRes) const = 0;

  virtual void
  computeForOneDomain(const int inDomainId, const util::rvec &box,
                      const energy::ForceField &forcefield,
                      const pairkernels::PairKernelManager &inKernels,
                      energy::EnergyMatrix<RealType> &outRes) const = 0;

  virtual std::string type() const = 0;
};

#endif
