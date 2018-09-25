// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#ifndef IO_CPLX_H
#define IO_CPLX_H
#include <map>
#include <yaml-cpp/yaml.h>

#include "domains/system.h"
#include "pairkernels/pairkernelmanager.h"

namespace io {

class YamlNode;

domains::System readCPLX(const std::string& file,
                         const std::vector<std::string>& beadTypes);
pairkernels::PairKernelManager readCPLXKernels(const std::string& yamlfile,
                                               const domains::Domains& domains);
util::rvec readBoxCPLX(const std::string& file);
energy::ForceField readForceFieldCPLX(const std::string& file);

// Internal functions made public for easy testing.
using ConnectionsMap = std::map<std::string, domains::Connections>;

domains::BeadChainID parseChainIDstr(const std::string& rawID);
std::vector<domains::BeadChainID> parseChainIds(const YamlNode& node);
std::vector<domains::Bead> parseBeads(
    const YamlNode& node, const std::vector<std::string>& beadTypes);
std::unique_ptr<domains::AbstractDomain> parseAbstractDomain(
    const YamlNode& definitionsDomainsNode, const YamlNode& node,
    const std::vector<std::string>& beadTypes,
    const ConnectionsMap& connections, const int id,
    const std::unordered_map<std::string, int>& typeIdsFromName);
ConnectionsMap parseConnectionsMap(const YamlNode& rawConnectionsNode,
                                   const YamlNode& rawDomainsNode);
}  // namespace io

#endif  // IO_CPLX_H
