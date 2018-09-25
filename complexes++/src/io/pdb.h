// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#ifndef IO_PDB_H
#define IO_PDB_H

#include <fstream>

#include "io/io.h"

namespace io {

void writePDB(std::basic_ostream<char>& out, const domains::Domains& model,
              const util::rvec& box, const int i,
              const std::vector<std::string>& beadTypes);

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
                     const util::rvec& box, const std::string& file);
  ~PDBReader();

  std::shared_ptr<domains::Domains> nextFrame() final;
  bool hasNextFrame() const final;

 private:
  void readPDBFrame();

  int m_frame;
  mutable std::ifstream m_pdb;
};
}  // namespace io

#endif  // IO_PDB_H
