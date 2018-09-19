// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#include <string>

#include "complexes_test.h"
#include "gtest/gtest.h"

#include "domains/abstractdomain.h"
#include "domains/rigiddomain.h"
#include "io/cplx.h"
#include "io/ffparameters.h"
#include "io/io.h"
#include "io/yaml.h"

TEST(CPLX_TEST, parseUnkownDomain) {
  const auto beadTypes = io::readBeadTypes(dataDir + "bead-types");
  EXPECT_THROW_MESSAGE(
      io::readTopologies(dataDir + "unkownDomain.cplx", beadTypes),
      std::invalid_argument, "foo --> Unkown Domain type, id 0\n");
}

TEST(CPLX_TEST, nonExistingFile) {
  const auto beadTypes = io::readBeadTypes(dataDir + "bead-types");
  EXPECT_THROW_MESSAGE(io::readTopologies("/tmp/foo.cplx", beadTypes),
                       std::invalid_argument,
                       "/tmp/foo.cplx <-- does not exist.\n");
}

TEST(CPLX_TEST, noConnection) {
  const auto emptyNode = io::YamlNode(YAML::Node());
  const auto con = io::parseConnectionsMap(emptyNode, emptyNode);
  EXPECT_EQ(0u, con.size());
}

TEST(CPLX_TEST, checkConnections) {
  const auto input = io::YamlNode(dataDir + "onlyConnections.yaml");
  const auto con =
      io::parseConnectionsMap(input["connections"], input["domains"]);
  EXPECT_EQ(2u, input["connections"].size());
  EXPECT_EQ(3u, input["domains"].size());
  EXPECT_EQ(input["domains"].size(), con.size());

  const auto dom1Connections = con.at("dom1");
  EXPECT_EQ(1u, dom1Connections.size());
  EXPECT_EQ(1, dom1Connections.at(0)->domainId());
  EXPECT_EQ(0, dom1Connections.at(0)->beadSelf());
  EXPECT_EQ(0, dom1Connections.at(0)->beadOther());

  const auto dom2Connections = con.at("dom2");
  EXPECT_EQ(2u, dom2Connections.size());
  EXPECT_EQ(0, dom2Connections.at(0)->domainId());
  EXPECT_EQ(0, dom2Connections.at(0)->beadSelf());
  EXPECT_EQ(0, dom2Connections.at(0)->beadOther());
  EXPECT_EQ(2, dom2Connections.at(1)->domainId());
  EXPECT_EQ(1, dom2Connections.at(1)->beadSelf());
  EXPECT_EQ(0, dom2Connections.at(1)->beadOther());

  const auto dom3Connections = con.at("dom3");
  EXPECT_EQ(1u, dom3Connections.size());
  EXPECT_EQ(1, dom3Connections.at(0)->domainId());
  EXPECT_EQ(0, dom3Connections.at(0)->beadSelf());
  EXPECT_EQ(1, dom3Connections.at(0)->beadOther());
}

TEST(CPLX_TEST, beadChainParsing) {
  const auto beadChain = io::parseChainIDstr("A 1");
  EXPECT_STREQ("A", beadChain.chain().c_str());
  EXPECT_EQ(1, beadChain.beadID());

  EXPECT_THROW_MESSAGE(io::parseChainIDstr("A"), std::invalid_argument,
                       "Invalid Bead Chain ID: A\n");
  EXPECT_THROW_MESSAGE(io::parseChainIDstr("A A"), std::invalid_argument,
                       "Invalid Bead Chain ID: A A\n");
  EXPECT_THROW_MESSAGE(io::parseChainIDstr("A 1 A"), std::invalid_argument,
                       "Invalid Bead Chain ID: A 1 A\n");

  const auto node = io::YamlNode(YAML::Load("[A 1, A 2, A 3]"));
  const auto chainIDs = io::parseChainIds(node);
  EXPECT_EQ(3u, chainIDs.size());
}

TEST(CPLX_TEST, parseRigid) {
  const auto beadTypes = io::readBeadTypes(dataDir + "bead-types");
  const auto node = io::YamlNode(dataDir + "rigidDomain.yaml");
  const auto domain = io::parseAbstractDomain(
      node["definitions"]["domains"], node["domains"], beadTypes,
      io::ConnectionsMap(), 0, {{"rigid", 0}});

  EXPECT_EQ(3, domain->nBeads());
  EXPECT_EQ(domains::findBeadID("ALA", beadTypes), domain->beads().at(0));
  EXPECT_EQ(domains::findBeadID("HIS", beadTypes), domain->beads().at(2));
  EXPECT_EQ(0, domain->charges().at(0));
  EXPECT_EQ(1, domain->charges().at(1));
  EXPECT_EQ(-1, domain->charges().at(2));
  EXPECT_STREQ("rigid", domain->name().c_str());

  auto const xyz = domain->xyz();
  for (auto i = 0; i < domain->nBeads(); ++i) {
    for (auto j = 0u; j < 3; ++j) {
      EXPECT_EQ(i, xyz(i, j));
    }
  }
}

TEST(CPLX_TEST, singleBeadDomain) {
  const auto beadTypes = io::readBeadTypes(dataDir + "bead-types");
  const auto node = io::YamlNode(dataDir + "singleBeadDomain.yaml");
  io::parseAbstractDomain(node["definitions"]["domains"], node["domains"],
                          beadTypes, io::ConnectionsMap(), 0, {{"rigid", 0}});
}

TEST(CPLX_TEST, singleDomain) {
  const auto beadTypes = io::readBeadTypes(dataDir + "bead-types");
  io::readCPLX(dataDir + "singleDomain.cplx", beadTypes);
}

TEST(CPLX_TEST, wrong_nBeads) {
  const auto beadTypes = io::readBeadTypes(dataDir + "bead-types");
  const auto node = io::YamlNode(dataDir + "wrong-nbeads.yaml");
  EXPECT_THROW_MESSAGE(
      io::parseAbstractDomain(node["definitions"]["domains"], node["domains"],
                              beadTypes, io::ConnectionsMap(), 0,
                              {{"rigid", 0}}),
      std::invalid_argument,
      "Number of beads in domain (id=0) does not match with vector size: "
      "charges.size() = 3, expected 4\n");
}

TEST(CPLX_TEST, domainIDOrder) {
  const auto beadTypes = io::readBeadTypes(dataDir + "bead-types");
  EXPECT_THROW_MESSAGE(
      io::readCPLX(dataDir + "wrong-domain-idx.cplx", beadTypes),
      std::invalid_argument,
      "Domain ids aren't increasing by single increments\n");
}

TEST(CPLX_TEST, wrongBeadConnectionIdx) {
  const auto beadTypes = io::readBeadTypes(dataDir + "bead-types");
  try {
    io::readCPLX(dataDir + "bad-connection-idx.cplx", beadTypes);
  } catch (std::invalid_argument& e) {
    const std::string expectedStr("Can't find bead '");
    fmt::print("{}\n", e.what());
    EXPECT_STREQ(expectedStr.c_str(),
                 std::string(e.what()).substr(0, expectedStr.length()).c_str());
  }
}

TEST(CPLX_TEST, sameMoveDiffClass) {
  const auto beadTypes = io::readBeadTypes(dataDir + "bead-types");
  EXPECT_NO_THROW(
      io::readTopologies(dataDir + "same-move-diff-class.cplx", beadTypes));
}

TEST(CPLX_TEST, sameMoveSameClassErr) {
  const auto beadTypes = io::readBeadTypes(dataDir + "bead-types");
  EXPECT_THROW_MESSAGE_CONTAINS(
      io::readTopologies(dataDir + "same-move-same-class-err.cplx", beadTypes),
      std::runtime_error, "No domainname 'rigiddom' found in definitions.\n");
}
