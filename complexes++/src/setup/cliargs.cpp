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
#include <iostream>
#include <string>

#include "parallelization/ompmanager.h"
#include "setup/cliargs.h"

namespace setup {


CLIArgs::CLIArgs(const int &argc, const char *const argv[]) : m_args() {
  CLI::App app{"COMPLEXES is a coarse grained simulation tool"};
  fillArgs(app);
  app.parse(argc,argv);
}

void CLIArgs::fillArgs(CLI::App& app){
    ADD_OPTION(std::vector<std::string>, multidir, {},
               "multiple simulation directories");
    ADD_OPTION(std::string, config, "", "config file");
    ADD_OPTION(bool, backup, true, "backup of files");
    ADD_OPTION(bool, rerun, false, "recalculate energies from trajectory");
    ADD_OPTION(std::string, restart, "", "restart file");
    ADD_OPTION(bool, version, false, "show version");
    ADD_OPTION(int, replex, 0, "number of sweeps between exchanges");
    ADD_OPTION(int, replex_stat, 1000,
               "number of sweeps between statistic output");
    ADD_OPTION(std::string, replex_accept, "", "exchange accept function");
    ADD_OPTION(std::string, movestats, "pertype",
               "specify the move statistics to show. Could be pertype, "
               "perdomain, all, none");
    ADD_OPTION(int, nb_threads, omp_get_max_threads(), "number of threads");
    ADD_OPTION(std::string, replex_verbosity, "stats",
               "exchange log verbosity (stats, all, none)");
}

std::string CLIArgs::value(const std::string &key) const {
  return value<std::string>(key);
}

bool CLIArgs::hasKey(const std::string &key) const noexcept {
  return static_cast<bool>(m_args.count(key));
}

} // namespace setup
