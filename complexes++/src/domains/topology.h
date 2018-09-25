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
#ifndef TOPOLOGY_H
#define TOPOLOGY_H

#include "domains/topologymoves.h"
#include "io/serializer.h"

#include <queue>
#include <set>
#include <unordered_map>
#include <vector>

namespace domains {
class Topology : public io::AbstractSerializable {
public:
  Topology(const std::vector<int> &inDomainIds, const bool inMove)
      : m_domainIds(inDomainIds), m_move(inMove), m_translationCoef(0),
        m_rotationCoef(0) {}

  const std::vector<int> &domainIds() const { return m_domainIds; }

  bool move(const util::rvec inBox, domains::Domains *outDomains,
            util::RNGEngine *inRng) {
    if (m_move == false) {
      return false;
    }

    auto moveIsEffective = false;
    auto dist = std::uniform_real_distribution<>{0, 1};
    if (dist(*inRng) < 0.5) {
      moveIsEffective = topologyMoves::Translation(
          m_domainIds, m_translationCoef, outDomains, inRng);
    } else {
      moveIsEffective = topologyMoves::Rotation(
          inBox, m_domainIds, m_rotationCoef, outDomains, inRng);
    }

    if (moveIsEffective) {
      for (auto domainId : m_domainIds) {
        util::pbc::applyPBCInPlace(inBox, &(*outDomains)[domainId]->xyz());
      }
    }

    return moveIsEffective;
  }

  void serialize(io::Serializer &serializer) const final {
    serializer.append(m_domainIds, "m_domainIds");
    serializer.append(m_move, "m_move");
  }

  Topology(io::Deserializer &deserializer)
      : m_domainIds(deserializer.restore<decltype(m_domainIds)>("m_domainIds")),
        m_move(deserializer.restore<decltype(m_move)>("m_move")) {}

  void setTranslationCoef(const double inTranslationCoef) {
    m_translationCoef = inTranslationCoef;
  }

  void setRotationCoef(const double inRotationCoef) {
    m_rotationCoef = inRotationCoef;
  }

  std::tuple<bool, int>
  allConnectionsAreValid(const domains::Domains &allDomains) const {
    if (m_domainIds.size() == 0) {
      return std::make_tuple(true, -1);
    }

    std::unordered_map<int, int> domainsIdsMapping;
    {
      int idxPosition = 0;
      for (const auto &dom : allDomains) {
        domainsIdsMapping[dom->id()] = idxPosition++;
      }
    }

    std::set<int> domainIdsSet;
    for (const auto id : m_domainIds) {
      if (domainsIdsMapping.find(id) == domainsIdsMapping.end()) {
        // this id does not exist
        return std::make_tuple(false, id);
      }
      domainIdsSet.insert(id);
    }

    for (auto id : m_domainIds) {
      const auto &connectionsForDomain =
          allDomains[domainsIdsMapping[id]]->connections();
      for (const auto &connection : connectionsForDomain) {
        const int connectionDomainId = connection->domainId();
        if (domainIdsSet.find(connectionDomainId) == domainIdsSet.end()) {
          // this connection is not part of the topology
          return std::make_tuple(false, id);
        }
      }
    }

    return std::make_tuple(true, -1);
  }

  std::tuple<bool, int>
  allDomainsConnected(const domains::Domains &allDomains) const {
    if (m_domainIds.size() == 0) {
      return std::make_tuple(true, -1);
    }

    std::unordered_map<int, int> domainsIdsMapping;
    {
      int idxPosition = 0;
      for (const auto &dom : allDomains) {
        domainsIdsMapping[dom->id()] = idxPosition++;
      }
    }

    std::set<int> domainIdsSet;
    for (const auto id : m_domainIds) {
      if (domainsIdsMapping.find(id) == domainsIdsMapping.end()) {
        // this id does not exist
        return std::make_tuple(false, id);
      }
      domainIdsSet.insert(id);
    }

    std::set<int> visitedDomainIds;
    std::queue<int> toVisitiDomainIds;
    toVisitiDomainIds.push(m_domainIds.front());

    while (!toVisitiDomainIds.empty()) {
      const auto currentDomainId = toVisitiDomainIds.front();
      toVisitiDomainIds.pop();

      if (visitedDomainIds.find(currentDomainId) == visitedDomainIds.end()) {
        visitedDomainIds.insert(currentDomainId);

        const auto &connectionsForDomain =
            allDomains[domainsIdsMapping[currentDomainId]]->connections();
        for (const auto &connection : connectionsForDomain) {
          const int connectionDomainId = connection->domainId();
          if (domainIdsSet.find(connectionDomainId) == domainIdsSet.end()) {
            // this connection id is not part of the topology
            return std::make_tuple(false, currentDomainId);
          }
          if (visitedDomainIds.find(connectionDomainId) ==
              visitedDomainIds.end()) {
            toVisitiDomainIds.push(connectionDomainId);
          }
        }
      }
    }

    if (visitedDomainIds == domainIdsSet) {
      return std::make_tuple(true, -1);
    } else {
      for (const auto id : m_domainIds) {
        if (visitedDomainIds.find(id) == visitedDomainIds.end()) {
          return std::make_tuple(false, id);
        }
      }
      return std::make_tuple(false, -1);
    }
  }

private:
  std::vector<int> m_domainIds;
  bool m_move;
  double m_translationCoef;
  double m_rotationCoef;
};
} // namespace domains
#endif // TOPOLOGY_H
