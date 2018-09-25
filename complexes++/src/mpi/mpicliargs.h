// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#ifndef MPI_MPICLIARGS_H
#define MPI_MPICLIARGS_H

#include <boost/program_options.hpp>
#include <fmt/format.h>

#include "setup/cliargs.h"

namespace mpi {
class MpiCLIArgs : public setup::CLIArgs {
 public:
  MpiCLIArgs(const int& argc, const char* const argv[]) : setup::CLIArgs() {
    addMpiArgs();
    boost::program_options::store(
        boost::program_options::command_line_parser(argc, argv)
            .options(m_required)
            .run(),
        m_args);
    boost::program_options::notify(m_args);
  }

 private:
  void addMpiArgs() {
    boost::program_options::options_description mpirequired("Requiered");

    mpirequired.add_options()(
        "mpi-partitions",
        boost::program_options::value<std::vector<int> >()->multitoken(),
        "number of simulation per nodes");

    m_required.add(mpirequired);
  }
};

}  // namespace config

#endif  // CONFIG_CLIARGS_H
