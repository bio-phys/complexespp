// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#ifndef MPI_MPIUTILS_H
#define MPI_MPIUTILS_H

#include <stdexcept>
// Do not use the openmpi c++ binding
#define OMPI_SKIP_MPICXX
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
}

#define MPI_ASSERT(COMMAND) mpi::mpiassert(COMMAND, __LINE__, __FILE__);

#endif
