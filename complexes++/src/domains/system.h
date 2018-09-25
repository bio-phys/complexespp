// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The coddomainIds(inDomainIds), m_move(inMove)domainIds(inDomainIds),
// m_move(inMove)e comes without warranty of any kind Please refer to Kim and
// Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#ifndef SYSTEM_H
#define SYSTEM_H

#include <memory>

#include "domains/abstractdomain.h"
#include "domains/topology.h"

namespace domains {
class System {
 public:
  System(std::shared_ptr<Domains> inDomains, std::vector<Topology> inTopologies)
      : m_domains(inDomains), m_topologies(inTopologies) {}

  std::shared_ptr<Domains> domains() { return m_domains; }
  const std::shared_ptr<Domains> domains() const { return m_domains; }
  std::vector<Topology> topologies() const { return m_topologies; }

 private:
  std::shared_ptr<Domains> m_domains;
  std::vector<Topology> m_topologies;
};
}  // namespace domains
#endif  // SYSTEM_H
