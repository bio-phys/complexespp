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
#ifndef IO_CPLX_H
#define IO_CPLX_H
#include <map>
#include <yaml-cpp/yaml.h>

#include "domains/system.h"
#include "pairkernels/pairkernelmanager.h"

namespace io {

class YamlNode;

domains::System readCPLX(const std::string &file,
                         const std::vector<std::string> &beadTypes);
pairkernels::PairKernelManager readCPLXKernels(const std::string &yamlfile,
                                               const domains::Domains &domains);
util::rvec readBoxCPLX(const std::string &file);
energy::ForceField readForceFieldCPLX(const std::string &file);

// Internal functions made public for easy testing.
using ConnectionsMap = std::map<std::string, domains::Connections>;

domains::BeadChainID parseChainIDstr(const std::string &rawID);
std::vector<domains::BeadChainID> parseChainIds(const YamlNode &node);
std::vector<domains::Bead>
parseBeads(const YamlNode &node, const std::vector<std::string> &beadTypes);
std::unique_ptr<domains::AbstractDomain> parseAbstractDomain(
    const YamlNode &definitionsDomainsNode, const YamlNode &node,
    const std::vector<std::string> &beadTypes,
    const ConnectionsMap &connections, const int id,
    const std::unordered_map<std::string, int> &typeIdsFromName);
ConnectionsMap parseConnectionsMap(const YamlNode &rawConnectionsNode,
                                   const YamlNode &rawDomainsNode);
} // namespace io

#endif // IO_CPLX_H
