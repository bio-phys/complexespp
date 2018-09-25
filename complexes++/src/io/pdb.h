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
#ifndef IO_PDB_H
#define IO_PDB_H

#include <fstream>

#include "io/io.h"

namespace io {

void writePDB(std::basic_ostream<char> &out, const domains::Domains &model,
              const util::rvec &box, const int i,
              const std::vector<std::string> &beadTypes);

// This implementation is not very smart. It doesn't do a lot of error
// checking.
// e.g. it won't notice any differences in the model and structures saved in
// the
// pdb as long as the number of atoms is correct. It pretty much only checks
// that the number of atoms is correct. Is doesn't check the atom type or if
// the
// chains correspond the to correct domain.
class PDBReader : public BaseReader {
public:
  explicit PDBReader(std::shared_ptr<domains::Domains> dom,
                     const util::rvec &box, const std::string &file);
  ~PDBReader();

  std::shared_ptr<domains::Domains> nextFrame() final;
  bool hasNextFrame() const final;

private:
  void readPDBFrame();

  int m_frame;
  mutable std::ifstream m_pdb;
};
} // namespace io

#endif // IO_PDB_H
