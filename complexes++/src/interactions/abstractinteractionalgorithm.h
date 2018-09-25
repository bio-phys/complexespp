// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
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
                               const std::shared_ptr<domains::Domains>&>,
      public io::AbstractSerializable {
 protected:
  void serializeCore(io::Serializer& serializer) const {
    serializer.append(type(), "type");
  }

 public:
  AbstractInteractionAlgorithm(io::Deserializer& deserializer) {}
  AbstractInteractionAlgorithm() {}
  virtual ~AbstractInteractionAlgorithm() {}

  virtual void updateDomain(const int inDomainId) = 0;

  virtual void resetAllDomains(const util::rvec& newBox) = 0;

  virtual void computeAll(const util::rvec& box,
                          const energy::ForceField& forcefield,
                          const pairkernels::PairKernelManager& inKernels,
                          energy::EnergyMatrix<RealType>& outRes) const = 0;

  virtual void computeForOneDomain(
      const int inDomainId, const util::rvec& box,
      const energy::ForceField& forcefield,
      const pairkernels::PairKernelManager& inKernels,
      energy::EnergyMatrix<RealType>& outRes) const = 0;

  virtual std::string type() const = 0;
};

#endif
