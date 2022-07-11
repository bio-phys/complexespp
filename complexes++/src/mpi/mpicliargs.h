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
#ifndef MPI_MPICLIARGS_H
#define MPI_MPICLIARGS_H
#define FMT_HEADER_ONLY

#include "setup/cliargs.h"

namespace mpi {
class MpiCLIArgs : public setup::CLIArgs {
public:
  MpiCLIArgs(const int &argc, const char *const argv[]) : setup::CLIArgs(argc, argv) {
      args.addMultiParameter<int>({"mpi-partitions", "np"}, "number of simulation per nodes",
                             mpi_partitions);
  }

  std::vector<int> mpi_partitions;
};

} // namespace mpi

#endif // CONFIG_CLIARGS_H
