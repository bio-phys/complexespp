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
#include <algorithm>
#include <fmt/format.h>
#include <fstream>
#include <iterator>
#include <utility>
#include <yaml-cpp/yaml.h>

// AUTO-ADD-DOMAIN HEADER (this comment must remain here to use "make
// add-domain"!)
#include "domains/membrane.h"
#include "domains/rigiddomain.h"
#include "io/cplx.h"
#include "io/yaml.h"
#include "membrane/abstractmembrane.h"
#include "membrane/flatmembrane.h"
#include "membrane/spheremembrane.h"
#include "membrane/tubemembrane.h"
#include "util/linalg.h"
#include "util/log.h"
#include "util/string.h"

namespace io {

util::rvec readBoxCPLX(const std::string &file) {
  return parse3DVector<double>(YamlNode(file)["box"]);
}

domains::BeadChainID parseChainIDstr(const std::string &rawID) {
  const auto tokens = util::splitStr(rawID, " ");
  if (tokens.size() != 2) {
    throw std::invalid_argument(
        fmt::format("Invalid Bead Chain ID: {}\n", rawID));
  }
  try {
    return domains::BeadChainID(tokens[0], std::stoi(tokens[1]));
  } catch (std::invalid_argument &e) {
    throw std::invalid_argument(
        fmt::format("Invalid Bead Chain ID: {}\n", rawID));
  }
}

std::vector<domains::BeadChainID> parseChainIds(const YamlNode &node) {
  const auto rawChainIDs = parseVector<std::string>(node);
  auto chainIDs = std::vector<domains::BeadChainID>();
  chainIDs.resize(rawChainIDs.size(), domains::BeadChainID("A", 0));
  std::transform(rawChainIDs.begin(), rawChainIDs.end(), chainIDs.begin(),
                 [](const auto &id) { return parseChainIDstr(id); });
  return chainIDs;
}

std::vector<domains::Bead>
parseBeads(const YamlNode &node, const std::vector<std::string> &beadTypes) {
  const auto rawBeads = parseVector<std::string>(node);
  auto beads = std::vector<domains::Bead>(rawBeads.size());
  std::transform(rawBeads.begin(), rawBeads.end(), beads.begin(),
                 [&beadTypes](const auto &b) {
                   return domains::findBeadID(b, beadTypes);
                 });
  return beads;
}

std::unique_ptr<domains::Rigid>
parseRigid(const std::string &inDomainTypename, const YamlNode &node,
           const std::vector<std::string> &beadTypes,
           const ConnectionsMap &connections, const int id,
           const YamlNode &parametersNode, const int domainTypeID) {
  const auto name = node["name"].as<std::string>();
  const auto phi = parametersNode["rotation"].as<double>();
  const auto trans = [](const auto &n) {
    switch (n.Type()) {
    case YAML::NodeType::Scalar:
      return util::rvec(n.template as<double>());
    default:
      return parse3DVector<double>(n);
    }
  }(parametersNode["translation"]);
  const auto beads = parseBeads(node["beads"], beadTypes);
  const auto charges = parseVector<double>(node["charges"]);
  const auto chainIDs = parseChainIds(node["chain-ids"]);
  auto con = domains::Connections(0);
  if (connections.find(name) != std::cend(connections)) {
    con = connections.at(name);
  }

  auto dom = std::make_unique<domains::Rigid>(inDomainTypename, domainTypeID,
                                              id, beads, charges, chainIDs, con,
                                              trans, phi);
  dom->setXyz(parseArray<double>(node["coordinates"]));
  return dom;
}

std::unique_ptr<domains::Membrane>
parseFlatMembrane(const std::string &inDomainTypename, const YamlNode &node,
                  const int id, const YamlNode &parametersNode,
                  const std::vector<std::string> &beadTypes,
                  const int domainTypeID) {
  const auto zaxis = node["zaxis"].as<double>();

  const auto z0 = parametersNode["z0"].as<double>();
  const auto ePSI0 = parametersNode["ePSI0"].as<double>();
  const auto beadID = domains::findBeadID("MMM", beadTypes);
  auto dom = std::make_unique<domains::Membrane>(
      inDomainTypename, std::make_unique<membrane::Flat<double>>(zaxis), z0,
      ePSI0, domainTypeID, id, beadID, domains::BeadChainID("Z", 1));
  return dom;
}

std::unique_ptr<domains::Membrane>
parseTubeMembrane(const std::string &inDomainTypename, const YamlNode &node,
                  const int id, const YamlNode &parametersNode,
                  const std::vector<std::string> &beadTypes,
                  const int domainTypeID) {
  const auto x = node["x"].as<double>();
  const auto y = node["y"].as<double>();
  const auto radius = node["radius"].as<double>();

  const auto z0 = parametersNode["z0"].as<double>();
  const auto ePSI0 = parametersNode["ePSI0"].as<double>();
  const auto beadID = domains::findBeadID("MMM", beadTypes);
  auto dom = std::make_unique<domains::Membrane>(
      inDomainTypename, std::make_unique<membrane::Tube<double>>(x, y, radius),
      z0, ePSI0, domainTypeID, id, beadID, domains::BeadChainID("Z", 1));
  return dom;
}

std::unique_ptr<domains::Membrane>
parseSphereMembrane(const std::string &inDomainTypename, const YamlNode &node,
                    const int id, const YamlNode &parametersNode,
                    const std::vector<std::string> &beadTypes,
                    const int domainTypeID) {
  const auto center = parse3DVector<double>(node["center"]);
  const auto radius = node["radius"].as<double>();

  const auto z0 = parametersNode["z0"].as<double>();
  const auto ePSI0 = parametersNode["ePSI0"].as<double>();
  const auto beadID = domains::findBeadID("MMM", beadTypes);
  auto dom = std::make_unique<domains::Membrane>(
      inDomainTypename,
      std::make_unique<membrane::Sphere<double>>(radius, center), z0, ePSI0,
      domainTypeID, id, beadID, domains::BeadChainID("Z", 1));
  return dom;
}

std::unique_ptr<domains::Membrane>
parseMembrane(const std::string &inDomainTypename, const YamlNode &node,
              const std::vector<std::string> &beadTypes,
              const ConnectionsMap &connections, const int id,
              const YamlNode &parametersNode, const int domainTypeID) {
  UNUSED(connections);
  const auto memType = parametersNode["type"].as<std::string>();

  if (memType == membrane::Flat<double>::Type()) {
    return parseFlatMembrane(inDomainTypename, node, id, parametersNode,
                             beadTypes, domainTypeID);
  } else if (memType == membrane::Tube<double>::Type()) {
    return parseTubeMembrane(inDomainTypename, node, id, parametersNode,
                             beadTypes, domainTypeID);
  } else if (memType == membrane::Sphere<double>::Type()) {
    return parseSphereMembrane(inDomainTypename, node, id, parametersNode,
                               beadTypes, domainTypeID);
  } else {
    throw std::invalid_argument(
        fmt::format("{} --> Unkown Membrane type, id {}\n", memType, id));
  }
}

std::unique_ptr<domains::AbstractDomain> parseAbstractDomain(
    const YamlNode &definitionsDomainsNode, const YamlNode &node,
    const std::vector<std::string> &beadTypes,
    const ConnectionsMap &connections, const int id,
    const std::unordered_map<std::string, int> &typeIdsFromName) {
  // AUTO-ADD-DOMAIN START-PARSE (this comment must remain here to use "make
  // add-domain"!)
  const auto typeNode = node["type"];

  const std::string domainTypename = typeNode.as<std::string>();
  auto domainTypeID = -1;
  try {
    domainTypeID = typeIdsFromName.at(domainTypename);
  } catch (const std::out_of_range &) {
    throw std::runtime_error(fmt::format(
        "No domainname '{}' found in definitions.\n", domainTypename));
  }

  const auto definitionForDom = definitionsDomainsNode[domainTypename];
  if (!definitionForDom.IsDefined()) {
    throw std::invalid_argument(fmt::format(
        "move of for definition {} cannot be found\n", domainTypename));
  }

  const auto definitionForDomMove = definitionForDom["move"];
  const std::string type = definitionForDomMove.as<std::string>();

  if (type == domains::Rigid::Type()) {
    return parseRigid(domainTypename, node, beadTypes, connections, id,
                      definitionForDom["defaults"], domainTypeID);
  } else if (type == domains::Membrane::Type()) {
    return parseMembrane(domainTypename, node, beadTypes, connections, id,
                         definitionForDom["defaults"], domainTypeID);
  }
  // AUTO-ADD-DOMAIN END-PARSE (this comment must remain here to use "make
  // add-domain"!)
  else {
    throw std::invalid_argument(
        fmt::format("{} --> Unkown Domain type, id {}\n", type, id));
  }
}

std::string formatNode(const YamlNode &node) {
  auto str = std::string("[");
  const auto size = static_cast<int>(node.size());
  for (auto i = 0; i < size - 1; ++i) {
    str += node[i].as<std::string>() + ", ";
  }
  str += node[size - 1].as<std::string>() + "]";
  return str;
}

std::tuple<std::string, int, int> parseConnectionDomain(
    const YamlNode &node, const std::map<std::string, int> &names2id,
    const std::map<std::string, std::vector<domains::BeadChainID>>
        &names2chainIds,
    const int conIdx) {
  const auto name = node[0].as<std::string>();
  auto id = -1;
  try {
    id = names2id.at(name);
  } catch (const std::out_of_range &e) {
    throw std::invalid_argument(fmt::format(
        "didn't found domain '{}' in topology when parsing connections\n",
        name));
  }
  const auto beadChainId = parseChainIDstr(node[1].as<std::string>());
  const auto beadId =
      util::indexOf(std::cbegin(names2chainIds.at(name)),
                    std::cend(names2chainIds.at(name)), beadChainId);
  if (beadId == -1) {
    throw std::invalid_argument(
        fmt::format("Can't find bead '{}' in domain {} in connection {}\n",
                    beadChainId, name, conIdx));
  }
  return std::make_tuple(name, id, beadId);
}

std::shared_ptr<domains::Connection>
parseConnection(const int bead, const int otherId, const int otherBead,
                const std::string &type, const YamlNode &node) {
  if (type == "flat") {
    return std::make_shared<domains::FlatConnection>(bead, otherId, otherBead);
  } else if (type == "harmonic") {
    return std::make_shared<domains::HarmonicConnection>(
        bead, otherId, otherBead, node["x0"].as<double>(),
        node["k"].as<double>());
  } else if (type == "gaussian") {
    return std::make_shared<domains::GaussianConnection>(
        bead, otherId, otherBead, node["N"].as<int>(),
        node["bond-length"].as<double>());
  }

  throw std::invalid_argument(
      fmt::format("Unknown connection type: {}\n", type));
}

ConnectionsMap parseConnectionsMap(const YamlNode &rawConnectionsNode,
                                   const YamlNode &rawDomainsNode) {
  if (!rawConnectionsNode.IsDefined() || rawConnectionsNode.size() == 0) {
    return ConnectionsMap();
  }

  if (rawConnectionsNode.Type() != YAML::NodeType::Map) {
    throw std::invalid_argument(
        fmt::format("Connections node is expected to be a map\n"));
  }

  auto namesToId = std::map<std::string, int>();
  auto namesToChainId =
      std::map<std::string, std::vector<domains::BeadChainID>>();
  const auto ids = rawDomainsNode.keys<int>();
  for (const auto i : ids) {
    const auto dom = rawDomainsNode[std::to_string(i)];
    const auto name = dom["name"].as<std::string>();
    namesToId[name] = i;
    namesToChainId[name] = parseChainIds(dom["chain-ids"]);
  }

  auto connections_map = ConnectionsMap();
  for (auto idx = 0u; idx < rawConnectionsNode.size(); ++idx) {
    const auto &con = rawConnectionsNode[std::to_string(idx)];

    try {
      const auto firstDom = parseConnectionDomain(con["domain-a"], namesToId,
                                                  namesToChainId, idx);
      const auto firstDomainName = std::get<0>(firstDom);
      const auto firstID = std::get<1>(firstDom);
      const auto firstBead = std::get<2>(firstDom);

      const auto secondDom = parseConnectionDomain(con["domain-b"], namesToId,
                                                   namesToChainId, idx);
      const auto secondDomainName = std::get<0>(secondDom);
      const auto secondID = std::get<1>(secondDom);
      const auto secondBead = std::get<2>(secondDom);

      const auto type = con["type"].as<std::string>();

      // insert first domain connection
      if (connections_map.find(firstDomainName) == std::end(connections_map)) {
        connections_map[firstDomainName] = domains::Connections();
      }
      connections_map[firstDomainName].emplace_back(parseConnection(
          firstBead, secondID, secondBead, type, con["params"]));

      // insert second domain connection
      if (connections_map.find(secondDomainName) == std::end(connections_map)) {
        connections_map[secondDomainName] = domains::Connections();
      }
      connections_map[secondDomainName].emplace_back(
          parseConnection(secondBead, firstID, firstBead, type, con["params"]));
    } catch (const YAML::InvalidNode &e) {
      throw e;
    }
  }

  return connections_map;
}

std::vector<int> parseDomainIDs(const YamlNode &node) {
  auto ids = node.keys<int>();

  // Sort the ids that have are only dictionary keys because I can't guarantee
  // the order in which they are read.
  std::sort(std::begin(ids), std::end(ids));
  if (!util::isConsecutive(ids)) {
    throw std::invalid_argument(
        "Domain ids aren't increasing by single increments\n");
  }

  return ids;
}

domains::Domains
parseDomains(const YamlNode &node, const YamlNode &definitionsDomainsNode,
             const std::vector<std::string> &beadTypes,
             const std::unordered_map<std::string, int> &typeIdsFromName) {
  const auto &domainNode = node["domains"];
  const auto nDomains = node["ndomains"].as<std::size_t>();
  if (nDomains == 0) {
    throw std::invalid_argument(fmt::format("No domain defined\n"));
  }
  if (domainNode.size() != nDomains) {
    throw std::invalid_argument(fmt::format("Specified {} domains, found: {}\n",
                                            nDomains, domainNode.size()));
  }
  const auto domainIDs = parseDomainIDs(domainNode);

  const auto connections = parseConnectionsMap(node["connections"], domainNode);

  auto topology = domains::Domains();
  for (const auto domID : domainIDs) {
    const auto dom = domainNode[domID];
    topology.emplace_back(parseAbstractDomain(definitionsDomainsNode, dom,
                                              beadTypes, connections, domID,
                                              typeIdsFromName));
  }

  return topology;
}

std::vector<std::array<std::string, 3>>
readPairInteraction(const YamlNode &pairInteractionNode) {
  if (pairInteractionNode.Type() != YAML::NodeType::Sequence) {
    throw std::invalid_argument(fmt::format(
        "Problem when reading pair-interaction in cplx, it is not a "
        "sequence (type {})\n",
        pairInteractionNode.Type()));
  }

  std::vector<std::array<std::string, 3>> allPairs;

  for (auto idx = 0u; idx < pairInteractionNode.size(); ++idx) {
    auto pair = pairInteractionNode[idx];

    if (pair.Type() != YAML::NodeType::Map) {
      throw std::invalid_argument(fmt::format(
          "Problem when reading pair-interaction in cplx, it is not a "
          "map at position {} (type {})\n",
          idx, pair.Type()));
    }

    const std::string kernelName = pair["function"].as<std::string>();

    auto pairNames = pair["domain-type-pair"];

    if (pairNames.Type() != YAML::NodeType::Sequence) {
      throw std::invalid_argument(
          fmt::format("Problem when reading pair-interaction in cplx, "
                      "domain-type-pair is not a "
                      "sequence at position {} (type {})\n",
                      idx, pairNames.Type()));
    }
    if (pairNames.size() != 2) {
      throw std::invalid_argument(fmt::format(
          "Problem when reading pair-interaction in cplx, domain-type-pair "
          "must contain 2 values, at position {} (size {})\n",
          idx, pairNames.size()));
    }

    const std::string domainType1 = pairNames[0].as<std::string>();
    const std::string domainType2 = pairNames[1].as<std::string>();

    allPairs.emplace_back(
        std::array<std::string, 3>{{domainType1, domainType2, kernelName}});
  }

  return allPairs;
}

domains::System readCPLX(const std::string &yamlfile,
                         const std::vector<std::string> &beadTypes) {
  TIMEZONE("read CPLX")
  const auto input = YamlNode(yamlfile);
  const auto &topologiesNode = input["topologies"];

  if (!topologiesNode.IsDefined()) {
    throw std::invalid_argument(fmt::format("No topology defined\n"));
  }

  const auto definitionsDomainsNode = input["definitions"]["domains"];
  std::unordered_map<std::string, int> typeIdsFromName;
  auto id = 0;
  for (const auto &name : definitionsDomainsNode.keys<std::string>()) {
    typeIdsFromName[name] = id++;
  }

  auto allDomains = domains::Domains();
  auto topologies = std::vector<domains::Topology>();
  // for (const auto topoNode : topologiesNode) {
  for (auto i = 0u; i < topologiesNode.size(); ++i) {
    const auto topologyNode = topologiesNode[i];
    auto top = parseDomains(topologyNode, definitionsDomainsNode, beadTypes,
                            typeIdsFromName);

    const auto full_mode = topologyNode["full-move"].as<bool>();

    // append to topology
    topologies.emplace_back(
        [](auto &t) {
          auto a = std::vector<int>();
          for (const auto &d : t) {
            a.push_back(d->id());
          }
          return a;
        }(top),
        full_mode);

    bool topologyIsValid;
    int badDomainId;
    std::tie(topologyIsValid, badDomainId) =
        topologies.back().allConnectionsAreValid(top);
    if (topologyIsValid == false) {
      throw std::invalid_argument(fmt::format(
          "All domains inside a topology must be connected with domains of the "
          "same topology. This is not the case for topology n° {} and domain "
          "id {}.\n",
          i, badDomainId));
    }

    bool topologyIsConnected;
    std::tie(topologyIsConnected, badDomainId) =
        topologies.back().allDomainsConnected(top);
    if (topologyIsConnected == false) {
      auto errorMessage =
          fmt::format("All domains inside a topology must be connected. This "
                      "is not the case for topology n° {} and domain id {}.\n",
                      i, badDomainId);
      if (full_mode) {
        throw std::invalid_argument(errorMessage);
      } else {
        std::cerr << errorMessage << std::endl;
      }
    }

    if (full_mode) {
      topologies.back().setTranslationCoef(
          topologyNode["translation"].as<double>());
      topologies.back().setRotationCoef(topologyNode["rotation"].as<double>());
    }

    allDomains.insert(allDomains.end(), std::make_move_iterator(top.begin()),
                      std::make_move_iterator(top.end()));
  }

  auto temp = std::make_shared<domains::Domains>(std::move(allDomains));
  return domains::System(temp, topologies);
}

pairkernels::PairKernelManager
readCPLXKernels(const std::string &yamlfile, const domains::Domains &domains) {
  TIMEZONE("read CPLX Kernels")
  const auto input = YamlNode(yamlfile);

  const auto &pairInteractionNode = input["definitions"]["pair-interaction"];
  auto allPairs = readPairInteraction(pairInteractionNode);

  return pairkernels::PairKernelManager(allPairs, domains);
}

energy::PairParameter<double>
readPairParameter(const YamlNode &ff, const std::vector<std::string> &beadTypes,
                  const std::string &name) {
  const auto node = ff[name];
  const auto nBeads = beadTypes.size();
  auto count = 0u;
  auto array = util::rArray(nBeads, nBeads);
  for (const auto &a : node.keys<std::string>()) {
    for (const auto &b : node[a].keys<std::string>()) {
      const auto value = node[a][b].as<double>();
      array(domains::findBeadID(a, beadTypes),
            domains::findBeadID(b, beadTypes)) = value;
      array(domains::findBeadID(b, beadTypes),
            domains::findBeadID(a, beadTypes)) = value;
      count++;
    }
  }
  if (count != nBeads * (nBeads - 1) / 2 + nBeads) {
    throw std::invalid_argument(
        fmt::format("{} --> Couldn't find values for all bead pairs.\n"
                    "Found {} and expected {}.",
                    name, count, nBeads * (nBeads - 1) / 2 + nBeads));
  }
  return energy::PairParameter<double>(array);
}

std::vector<std::array<double, 8>>
readMembranePotential(const YamlNode &ff,
                      const std::vector<std::string> &beadTypes) {
  const auto node = ff["membrane"];
  const auto nBeads = beadTypes.size();
  if (node.size() != nBeads) {
    throw std::invalid_argument(
        fmt::format("membrane <-- couldn't find values for all beads.\n"
                    "Found {} and expected {}.",
                    node.size(), nBeads));
  }

  auto mem = std::vector<std::array<double, 8>>(nBeads);
  for (const auto &beadType : node.keys<std::string>()) {
    const auto a = parseVector<double>(node[beadType]);
    std::copy(a.begin(), a.end(),
              mem[domains::findBeadID(beadType, beadTypes)].begin());
  }
  return mem;
}

std::vector<double> readChargeRadii(const YamlNode &ff,
                                    const std::vector<std::string> &beadTypes) {
  const auto node = ff["charge-radii"];
  const auto nBeads = beadTypes.size();
  if (node.size() != nBeads) {
    throw std::invalid_argument(
        fmt::format("chargeRadii <-- couldn't find values for all beads.\n"
                    "Found {} and expected {}.",
                    node.size(), nBeads));
  }
  auto radii = std::vector<double>(nBeads);
  for (const auto &beadType : node.keys<std::string>()) {
    radii[domains::findBeadID(beadType, beadTypes)] =
        node[beadType].as<double>();
  }
  return radii;
}

energy::ForceField readForceFieldCPLX(const std::string &yamlfile) {
  const auto ff = YamlNode(yamlfile)["forcefield"];
  const auto beadTypes = parseVector<std::string>(ff["bead-types"]);
  const auto interActionEnergy = readPairParameter(ff, beadTypes, "energies");
  const auto diameter = readPairParameter(ff, beadTypes, "diameter");
  const auto membrane = readMembranePotential(ff, beadTypes);
  const auto chargeRadius = readChargeRadii(ff, beadTypes);
  return energy::ForceField(
      beadTypes, interActionEnergy, diameter, chargeRadius, membrane,
      ff["debye-length"].as<double>(), ff["dielectric-constant"].as<double>(),
      ff["alpha"].as<double>());
}

} // namespace io
