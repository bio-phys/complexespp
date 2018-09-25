// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#include <fmt/format.h>
#include <fstream>
#include <sstream>

#include "domains/beads.h"
#include "io/ffparameters.h"
#include "util/util.h"

namespace io {

energy::PairParameter<double> readPairParameter(
    const std::string& file, const std::vector<std::string>& beadTypes) {
  const auto nResidues = beadTypes.size();
  auto array = util::rArray(nResidues, nResidues);
  auto lineCount = 0u;
  std::ifstream input(file);

  std::string a, b, line;
  double value;

  while (std::getline(input, line)) {
    std::istringstream iss(line);
    iss >> a >> b >> value;
    array(domains::findBeadID(a, beadTypes),
          domains::findBeadID(b, beadTypes)) = value;
    array(domains::findBeadID(b, beadTypes),
          domains::findBeadID(a, beadTypes)) = value;
    ++lineCount;
  }

  // check that we read the whole upper trianlge matrix plus trace.
  if (lineCount != nResidues * (nResidues - 1) / 2 + nResidues) {
    throw std::invalid_argument(fmt::format(
        "{} <-- Couldn't find values for all bead pairs.\n"
        "Found {} and expected {}.",
        file, lineCount, nResidues * (nResidues - 1) / 2 + nResidues));
  }

  return energy::PairParameter<double>(array);
}

std::vector<std::array<double, 8>> readMembranePotential(
    const std::string& file, const std::vector<std::string>& beadTypes) {
  const auto nBeads = beadTypes.size();
  auto mem = std::vector<std::array<double, 8>>(nBeads);
  auto lineCount = 0u;
  std::ifstream input(file);
  std::string name, line;

  while (std::getline(input, line)) {
    std::istringstream iss(line);
    std::array<double, 8> a;
    iss >> name >> a[0] >> a[1] >> a[2] >> a[3] >> a[4] >> a[5] >> a[6] >> a[7];
    mem[domains::findBeadID(name, beadTypes)] = a;
    ++lineCount;
  }
  if (lineCount != nBeads) {
    throw std::invalid_argument(
        fmt::format("{} <-- couldn't find values for all beads.\n"
                    "Found {} and expected {}.",
                    file, lineCount, nBeads));
  }
  return mem;
}

std::vector<std::string> readBeadTypes(const std::string& file) {
  std::ifstream input(file);
  auto beadTypes = std::vector<std::string>();
  std::string line, beadCode;
  while (std::getline(input, line)) {
    std::istringstream iss(line);
    iss >> beadCode;
    if (beadCode.size() > 3 && beadCode.size() == 0) {
      throw std::invalid_argument(fmt::format(
          "Residue Code should be between 1 and 3 charaters. You used '{}'\n",
          beadCode));
    }
    beadTypes.push_back(beadCode);
  }
  return beadTypes;
}
}  // namespace io
