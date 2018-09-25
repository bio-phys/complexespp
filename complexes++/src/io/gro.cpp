// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#include <algorithm>
#include <fmt/format.h>
#include <fstream>

#include "domains/beads.h"
#include "io/gro.h"
#include "util/string.h"

namespace io {

namespace gro {
void writeATOMLine(std::basic_ostream<char>& out, const double x,
                   const double y, const double z, const std::string& resName,
                   const int atomID, const int resid) {
  const auto groFormat = "{:>5d}{:<5.5s}{:>5.5s}{:>5d}{:8.3f}{:8.3f}{:8.3f}\n";
  out << fmt::format(groFormat, util::truncateLeft(resid), resName, "CA",
                     util::truncateLeft(atomID), x, y, z);
}

void writeDom(const std::unique_ptr<domains::AbstractDomain>& dom,
              std::basic_ostream<char>& out,
              const std::vector<std::string>& beadTypes,
              const int atom_offset) {
  const auto nAtoms = dom->nBeads();
  for (auto i = 0; i < nAtoms; ++i) {
    const auto x = dom->xyz()(i, 0);
    const auto y = dom->xyz()(i, 1);
    const auto z = dom->xyz()(i, 2);
    const auto resName = beadTypes.at(dom->beads().at(i));
    const auto chain = dom->BeadChainIDs().at(i).chain();
    const auto resid = dom->BeadChainIDs().at(i).beadID();
    writeATOMLine(out, x, y, z, resName, i + atom_offset, resid);
  }
}
}  // namespace gro

void writeGRO(const std::string& fname, const domains::Domains& model,
              const util::rvec& box,
              const std::vector<std::string>& beadTypes) {
  std::ofstream out(fname);
  out << fmt::format("COMPLEXES STRUCTURE \n");
  const auto nbeads = std::accumulate(
      std::begin(model), std::end(model), 0,
      [](auto sum, const auto& dom) { return sum + dom->nBeads(); });
  out << fmt::format("{:>d}\n", nbeads);

  auto processed_beads = 0;
  for (auto const& dom : model) {
    gro::writeDom(dom, out, beadTypes, processed_beads);
    processed_beads += dom->nBeads();
  }
  out << fmt::format("{:10.5f}{:10.5f}{:10.5f}\n", box[0], box[1], box[2]);
}

}  // namespace io
