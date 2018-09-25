// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#include "domains/abstractdomain.h"
#include "domains/rigiddomain.h"

#include "complexes_test.h"
#include "gtest/gtest.h"

#include "io/pdb.h"

std::unique_ptr<domains::Rigid> dummyRidgidDomain(const int id,
                                                  const int beadID) {
  const auto nbeads = 4;
  const auto Dtrans = util::rvec(.2, .2, .2);
  const auto phi = 2.0;

  const auto beadChainIDs =
      std::vector<domains::BeadChainID>(nbeads, domains::BeadChainID("A", 0));

  auto dom = std::make_unique<domains::Rigid>(
      "none", 0, id, std::vector<domains::Bead>(nbeads, beadID),
      std::vector<double>(nbeads, 0), beadChainIDs, domains::Connections(0),
      Dtrans, phi);
  return dom;
}

domains::Domains dummyDomains(const int ndomains) {
  auto topology = domains::Domains();
  for (auto i = 0; i < ndomains; ++i) {
    topology.emplace_back(dummyRidgidDomain(i, i));
  }
  return topology;
}

TEST(PDB, writing) {
  // check if a PDB is written correct
  // I can read again from this after I have written to it
  std::stringstream to_write_stream;

  // create dummy domains
  const auto beadTypes = std::vector<std::string>{"ALA", "VAL"};
  const auto topology = dummyDomains(static_cast<int>(beadTypes.size()));
  const auto box = util::rvec(10, 10, 10);
  const auto frame = 0;

  io::writePDB(to_write_stream, topology, box, frame, beadTypes);

  // inspect written pdb
  std::string line;

  std::vector<std::string> expected_line = {
      {"CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1           1",
       "MODEL 1",
       "ATOM      1 CA   ALA A   0       0.000   0.000   0.000  0.00  0.00",
       "ATOM      2 CA   ALA A   0       0.000   0.000   0.000  0.00  0.00",
       "ATOM      3 CA   ALA A   0       0.000   0.000   0.000  0.00  0.00",
       "ATOM      4 CA   ALA A   0       0.000   0.000   0.000  0.00  0.00",
       "ATOM      5 CA   VAL A   0       0.000   0.000   0.000  0.00  0.00",
       "ATOM      6 CA   VAL A   0       0.000   0.000   0.000  0.00  0.00",
       "ATOM      7 CA   VAL A   0       0.000   0.000   0.000  0.00  0.00",
       "ATOM      8 CA   VAL A   0       0.000   0.000   0.000  0.00  0.00",
       "TER", "ENDMDL"}};

  auto i = 0;
  while (std::getline(to_write_stream, line)) {
    EXPECT_STREQ(expected_line.at(i).c_str(), line.c_str());
    ++i;
  }
}

TEST(PDB, very_large_system) {
  // check if a PDB is written correct
  // I can read again from this after I have written to it
  std::stringstream to_write_stream;

  // create dummy domains
  const auto beadTypes = std::vector<std::string>{"ALA", "VAL"};
  auto topology = domains::Domains();
  // a dummy domain has 4 beads. This will create 110004 beads in total
  for (auto i = 0; i < (110004 / 4); ++i) {
    topology.emplace_back(dummyRidgidDomain(i, 0));
  }
  const auto box = util::rvec(10, 10, 10);
  const auto frame = 0;

  io::writePDB(to_write_stream, topology, box, frame, beadTypes);

  std::string lastAtomLine;
  std::string line;

  while (std::getline(to_write_stream, line)) {
    if (line.substr(0, 4) == "ATOM") {
      lastAtomLine = line;
    }
  }
  // The first digit of the bead number will be truncated for ids above 99999.
  EXPECT_STREQ(
      "ATOM  10004 CA   ALA A   0       0.000   0.000   0.000  0.00  0.00",
      lastAtomLine.c_str());
}
