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
#define FMT_HEADER_ONLY
#include <fmt/format.h>
#include <sstream>

#include "domains/beads.h"
#include "io/pdb.h"
#include "util/string.h"

namespace io {

void writePDBLine(std::basic_ostream<char> &out, const double x, const double y,
                  const double z, const std::string &resName, const int atomID,
                  const std::string &chain, const int resid) {
  // taken from http://cupnet.net/pdb-format/
  const auto pdbFormat = "{:6s}{:5d} {:4s}{:1s}{:3s} {}{:4d}{:1s}   "
                         "{:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}\n";
  // shift atom ID by one because PDB is one indexed
  out << fmt::format(pdbFormat, "ATOM", util::truncateLeft(atomID + 1), "CA",
                     "", resName, chain, resid, "", x, y, z, 0.0, 0.0);
}

void writeDom(const std::unique_ptr<domains::AbstractDomain> &dom,
              std::basic_ostream<char> &out,
              const std::vector<std::string> &beadTypes,
              const int atom_offset) {
  const auto nAtoms = dom->nBeads();
  for (auto i = 0; i < nAtoms; ++i) {
    const auto x = dom->xyz()(i, 0);
    const auto y = dom->xyz()(i, 1);
    const auto z = dom->xyz()(i, 2);
    const auto resName = beadTypes.at(dom->beads().at(i));
    const auto chain = dom->BeadChainIDs().at(i).chain();
    const auto resid = dom->BeadChainIDs().at(i).beadID();
    writePDBLine(out, x, y, z, resName, i + atom_offset, chain, resid);
  }
}

void writeBox(std::basic_ostream<char> &out, const util::rvec &box) {
  const auto boxFormat =
      "CRYST1{:9.3f}{:9.3f}{:9.3f}{:7.2f}{:7.2f}{:7.2f} {:14s}1\n";
  const auto angle = 90.0;
  // The P 1 says that we are using a primitive box
  out << fmt::format(boxFormat, box[0], box[1], box[2], angle, angle, angle,
                     "P 1");
}

void writePDB(std::basic_ostream<char> &out, const domains::Domains &model,
              const util::rvec &box, const int i,
              const std::vector<std::string> &beadTypes) {
  writeBox(out, box);
  out << fmt::format("MODEL {}\n", i + 1);
  auto processed_beads = 0;
  for (auto const &dom : model) {
    writeDom(dom, out, beadTypes, processed_beads);
    processed_beads += dom->nBeads();
  }
  out << fmt::format("TER\nENDMDL\n");
}

PDBReader::PDBReader(std::shared_ptr<domains::Domains> dom,
                     const util::rvec &box, const std::string &file)
    : BaseReader(dom, box), m_frame(0), m_pdb(file) {}

PDBReader::~PDBReader() {}

util::rvec readAtomLine(const std::string &line) {
  return util::rvec(std::stod(line.substr(31, 7)),
                    std::stod(line.substr(39, 7)),
                    std::stod(line.substr(47, 7)));
}

util::rvec readNextAtom(std::ifstream &instream) {
  std::string line;
  while (std::getline(instream, line)) {
    if (line.substr(0, 4) == "ATOM") {
      return readAtomLine(line);
    } else if (line.substr(0, 6) == "ENDMDL") {
      throw std::invalid_argument(fmt::format(
          "Discoverd end of PDB model unexpected. The trajectory "
          "reader was still trying to populate the model with a new "
          "frame. Check that the model fits to the trajectory in the "
          "config file.\n"));
    }
  }
  throw std::invalid_argument(fmt::format(
      "Didn't find any ATOM in PDB Model. Check that the model fits to the "
      "trajectory in the config file.\n"));
}

void PDBReader::readPDBFrame() {
  for (auto &dom : *m_dom) {
    auto xyz = dom->xyz();
    for (auto i = 0; i < dom->nBeads(); ++i) {
      auto atom = readNextAtom(m_pdb);
      xyz(i, 0) = atom[0];
      xyz(i, 1) = atom[1];
      xyz(i, 2) = atom[2];
    }
    dom->setXyz(std::move(xyz));
  }

  auto bReachedEndOfModel = false;
  try {
    readNextAtom(m_pdb);
  } catch (std::invalid_argument &e) {
    bReachedEndOfModel = true;
  }

  if (!bReachedEndOfModel) {
    throw std::invalid_argument(fmt::format(
        "Could still read an atom after repopulating the model. Check that the "
        "model fits with the trajectory in the config file.\n"));
  }
}

bool PDBReader::hasNextFrame() const {
  // Checks if we find another 'MODEL' record
  const auto curPos = m_pdb.tellg();
  std::string line;
  while (std::getline(m_pdb, line)) {
    if (line.substr(0, 5) == "MODEL") {
      m_pdb.seekg(curPos);
      return true;
    }
  }
  return false;
}

std::shared_ptr<domains::Domains> PDBReader::nextFrame() {
  ++m_frame;
  readPDBFrame();
  return m_dom;
}
} // namespace io
