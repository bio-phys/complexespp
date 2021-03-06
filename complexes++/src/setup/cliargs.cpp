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

namespace po = boost::program_options;

namespace setup {

void printHelp(const CLIArgs &args) {
  std::cout << "COMPLEXES is a corse grained simulation tool\n\n";
  std::cout << "Description of CLI arguments" << std::endl;
  std::cout << std::endl;
  std::cout << args << std::endl;
}

CLIArgs::CLIArgs(const int &argc, const char *const argv[])
    : m_required("Required"), m_args() {
  defineArgs();
  po::store(po::command_line_parser(argc, argv).options(m_required).run(),
            m_args);
  po::notify(m_args);
}

CLIArgs::CLIArgs() : m_required("Required"), m_args() { defineArgs(); }

void CLIArgs::defineArgs() {
  m_required.add_options()("help,h", "print this help")(
      "config,c", po::value<std::string>(), "config file")(
      "multidir", po::value<std::vector<std::string>>()->multitoken(),
      "multi-simulation source directories")(
      "backup", po::value<bool>()->default_value(true),
      "turn on creation of backup files")(
      "rerun", po::value<bool>()->default_value(false),
      "recalculate energies from trajectory")(
      "restart", po::value<std::string>(), "restart file")(
      "version", "show version")("replex", po::value<int>(),
                                 "number of sweeps between exchanges")(
      "replex-stat", po::value<int>()->default_value(1000),
      "number of sweeps between statistic output")(
      "replex-accept", po::value<std::string>(), "exchange accept function")(
      "movestats", po::value<std::string>()->default_value("pertype"),
      "specify the move statistics to show. Could be pertype, perdomain, all, "
      "none")("nb-threads",
              po::value<int>()->default_value(omp_get_max_threads()),
              "number of threads")(
      "replex-verbosity", po::value<std::string>()->default_value("stats"),
      "exchange log verbosity (stats, all, none)");
}

std::string CLIArgs::value(const std::string &key) const {
  return value<std::string>(key);
}

bool CLIArgs::hasKey(const std::string &key) const noexcept {
  return static_cast<bool>(m_args.count(key));
}

std::ostream &CLIArgs::print(std::ostream &out) const {
  return out << m_required;
}

std::ostream &operator<<(std::ostream &out, const CLIArgs &args) {
  return args.print(out);
}

} // namespace setup
