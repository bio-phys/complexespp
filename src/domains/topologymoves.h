// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#ifndef TOPOLOGYMOVES_H
#define TOPOLOGYMOVES_H

#include "domains/abstractdomain.h"
#include "domains/topology.h"
#include "util/moves.h"
#include "util/pbc.h"
#include "util/quaternions/quat.h"

#include <queue>
#include <set>
#include <unordered_map>
#include <vector>

namespace domains {

namespace topologyMoves {

inline bool Translation(const std::vector<int>& inDomainIds,
                        const double inTranslationCoef,
                        domains::Domains* outDomains, util::RNGEngine* inRng) {
  if (inTranslationCoef == 0 || inDomainIds.size() == 0) {
    return false;
  }

  const double translationToApply =
      std::uniform_real_distribution<>{-1, 1}(*inRng) * inTranslationCoef;

  for (auto domainId : inDomainIds) {
    util::translateByConstantInPlace(translationToApply,
                                     &(*outDomains)[domainId]->xyz());
  }

  return true;
}

inline util::rvec CompactTopology(const util::rvec inSimulationBoxSize,
                                  const std::vector<int>& inDomainIds,
                                  domains::Domains* outDomains) {
  const long int nbDomainsInTopology =
      static_cast<long int>(inDomainIds.size());
  if (nbDomainsInTopology == 1) {
    return util::centroid((*outDomains)[inDomainIds[0]]->xyz());
  }

  struct DomainInfo {
    util::rvec centroid;
    util::rvec shift;
    int globalId;
  };

  std::vector<DomainInfo> domainInfos(nbDomainsInTopology);
  std::unordered_map<int, int> globalToLocalId;

  {
    long int idxIds = 0;
    for (auto domainId : inDomainIds) {
      domainInfos[idxIds].centroid =
          util::centroid((*outDomains)[domainId]->xyz());
      domainInfos[idxIds].globalId = domainId;
      domainInfos[idxIds].shift = util::rvec(0, 0, 0);
      globalToLocalId[domainId] = idxIds;

      idxIds += 1;
    }
  }

  std::set<int> visitedDomainIds;
  std::queue<int> toProcessDomainIds;

  toProcessDomainIds.push(inDomainIds[0]);

  while (!toProcessDomainIds.empty()) {
    const auto currentDomainGlobalId = toProcessDomainIds.front();
    const auto currentDomainLocalId = globalToLocalId[currentDomainGlobalId];
    toProcessDomainIds.pop();

    const auto& connectionsForDomain =
        (*outDomains)[currentDomainGlobalId]->connections();
    for (const auto& connection : connectionsForDomain) {
      const auto connectionDomainGlobalId = connection->domainId();
      if (visitedDomainIds.find(connectionDomainGlobalId) ==
          visitedDomainIds.end()) {
        visitedDomainIds.insert(connectionDomainGlobalId);
        toProcessDomainIds.push(connectionDomainGlobalId);

        const auto connectionDomainLocalId =
            globalToLocalId[connectionDomainGlobalId];
        util::rvec shiftDomTest;
        util::pbc::DistSquareBetweenPoints(
            domainInfos[currentDomainLocalId].centroid,
            domainInfos[connectionDomainLocalId].centroid, inSimulationBoxSize,
            &shiftDomTest);

        domainInfos[connectionDomainLocalId].shift = shiftDomTest;
        domainInfos[connectionDomainLocalId].centroid += shiftDomTest;
      }
    }
  }

  util::rvec currentCentroid(0, 0, 0);

  for (const auto& domainInfo : domainInfos) {
    util::translateByConstantInPlace(
        domainInfo.shift, &(*outDomains)[domainInfo.globalId]->xyz());
    currentCentroid += domainInfo.centroid;
  }

  for (int idxDim = 0; idxDim < 3; ++idxDim) {
    currentCentroid[idxDim] /= static_cast<double>(nbDomainsInTopology);
  }

  return currentCentroid;
}

inline bool Rotation(const util::rvec inSimulationBoxSize,
                     const std::vector<int>& inDomainIds,
                     const double inRotationCoef, domains::Domains* outDomains,
                     util::RNGEngine* inRng) {
  if (inRotationCoef == 0 || inDomainIds.size() == 0) {
    return false;
  }

  const auto rot = util::randomRotation<double>(inRotationCoef, *inRng);

  const auto com =
      CompactTopology(inSimulationBoxSize, inDomainIds, outDomains);

  for (auto domainId : inDomainIds) {
    util::rotateWithMatInPlace(com, rot, &(*outDomains)[domainId]->xyz());
  }

  return true;
}

}  // namespace topologyMoves

}  // namespace domains
#endif  // TOPOLOGY_H
