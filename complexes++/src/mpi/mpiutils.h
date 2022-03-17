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
#ifndef MPI_MPIUTILS_H
#define MPI_MPIUTILS_H

#include <stdexcept>
// Do not use the openmpi c++ binding
#define OMPI_SKIP_MPICXX
#define FMT_HEADER_ONLY
#include <fmt/format.h>
#include <mpi.h>

namespace mpi {

inline void mpiassert(const int mpiCommandRes, const int line,
                      const std::string file) {
  if (mpiCommandRes != MPI_SUCCESS) {
    throw std::runtime_error(fmt::format(
        "Error during the mpi call at line {} in file {}.", line, file));
  }
}
} // namespace mpi

#define MPI_ASSERT(COMMAND) mpi::mpiassert(COMMAND, __LINE__, __FILE__);

#endif
