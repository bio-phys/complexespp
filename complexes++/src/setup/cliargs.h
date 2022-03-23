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
#ifndef SETUP_CLIARGS_H
#define SETUP_CLIARGS_H

#define FMT_HEADER_ONLY
#include <clsimple.hpp>
#include <any>
#include <fmt/format.h>
#include <initializer_list>
#include <unordered_map>
#include <utility>

#include "parallelization/ompmanager.h"

namespace setup {
class CLIArgs {
public:
    CLIArgs(const int &argc, const char *const argv[]):
        args("COMPLEXES is a coarse grained simulation tool",
             argc, argv){

        args.addMultiParameter<std::string>("multidir", "multiple simulation directories",
                                        multidir);
        args.addParameter<std::string>("config", "config file",
                                       config, "");
        args.addParameter<bool>("backup", "backup of files",
                                backup, true);
        args.addParameter<bool>("rerun", "recalculate energies from trajectory",
                                rerun, false);
        args.addParameter<std::string>("restart", "restart file",
                                       restart, "");
        args.addParameter<bool>("version", "show version",
                                version, false);
        args.addParameter<int>("replex", "number of sweeps between exchanges",
                               replex, 0);
        args.addParameter<int>("replex_stat", "number of sweeps between statistic output",
                               replex_stat, 1000);
        args.addParameter<std::string>("replex_accept", "exchange accept function",
                                       replex_accept, "");
        args.addParameter<std::string>("movestats",
                                     "specify the move statistics to show. Could be pertype, "
                                     "perdomain, all, none",
                                      movestats, "pertype");
        args.addParameter<int>("nb_threads", "number of threads",
                               nb_threads, omp_get_max_threads());
        args.addParameter<std::string>("replex_verbosity", "exchange log verbosity (stats, all, none)",
                                       replex_verbosity, "stats");

    }

    bool hasKey(const std::string &key) const {
        return args.hasKey(key);
    }

    bool parse(){
      return args.parse();
    }

    template <class StreamClass>
    void printHelp(StreamClass& inStream) const{
        args.printHelp(inStream);
    }

    template <class ValType>
    static ValType getMapping(const std::string& inKey,
                              std::initializer_list<std::pair<std::string, ValType>> inMapping,
                              const ValType defaultVal = ValType()){
        return CLsimple::GetMapping(inKey, inMapping, defaultVal);
    }

    std::vector<std::string> multidir;
    std::string config;
    bool backup;
    bool rerun;
    std::string restart;
    bool version;
    int replex;
    int replex_stat;
    std::string replex_accept;
    std::string movestats;
    int nb_threads;
    std::string replex_verbosity;
protected:
  CLsimple args;

};

} // namespace setup

#endif // SETUP_CLIARGS_H
