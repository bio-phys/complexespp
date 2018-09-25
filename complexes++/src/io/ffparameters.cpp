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
#include <fmt/format.h>
#include <fstream>
#include <sstream>

#include "domains/beads.h"
#include "io/ffparameters.h"
#include "util/util.h"

namespace io {

energy::PairParameter<double>
readPairParameter(const std::string &file,
                  const std::vector<std::string> &beadTypes) {
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

std::vector<std::array<double, 8>>
readMembranePotential(const std::string &file,
                      const std::vector<std::string> &beadTypes) {
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

std::vector<std::string> readBeadTypes(const std::string &file) {
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
} // namespace io
