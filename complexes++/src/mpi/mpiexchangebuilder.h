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
#ifndef MPIEXCHANGEBUILDER_H
#define MPIEXCHANGEBUILDER_H

#include <memory>
#include <string>
#include <vector>

#include "mc/abstractexchangesimulation.h"
#include "mc/exchangeaccepter.h"
#include "mc/exchangelooper.h"
#include "mc/simulationchecking.h"
#include "mpi/mpibasicexchangesimulation.h"

namespace mpi {

std::unique_ptr<mc::AbstractExchangeSimulation>
MpiExchangeBuilder(const std::vector<std::string> &configDirNames,
                   const std::vector<std::string> &configDirNamesGlobal,
                   const std::string &configFilename, const int inExchangeRate,
                   const int inStatsRate, const std::string &accepterChoice,
                   const bool restart, const std::string &restartValue,
                   const bool backupOutput,
                   const io::MoveStatRecorder::Verbosity inMoveStatsVerbosity,
                   const int nbThreads, const std::vector<int> &partitions,
                   const std::vector<int> &partitionsOffset) {
  if (accepterChoice == "hrex") {
    return std::unique_ptr<mc::AbstractExchangeSimulation>(
        new MpiBasicExchangeSimulation<mc::OddEvenNeighborLoop, mc::HREXAccept>(
            configDirNames, configDirNamesGlobal, configFilename,
            inExchangeRate, inStatsRate, mc::SimulationChecking::ExchangeMask(),
            true, restart, restartValue, backupOutput, inMoveStatsVerbosity,
            nbThreads, partitions, partitionsOffset));
  } else if (accepterChoice == "remc") {
    return std::unique_ptr<mc::AbstractExchangeSimulation>(
        new MpiBasicExchangeSimulation<mc::OddEvenNeighborLoop, mc::REMCAccept>(
            configDirNames, configDirNamesGlobal, configFilename,
            inExchangeRate, inStatsRate, mc::SimulationChecking::ExchangeMask(),
            true, restart, restartValue, backupOutput, inMoveStatsVerbosity,
            nbThreads, partitions, partitionsOffset));
  } else if (accepterChoice == "npt") {
    return std::unique_ptr<mc::AbstractExchangeSimulation>(
        new MpiBasicExchangeSimulation<mc::OddEvenNeighborLoop, mc::NPTAccept>(
            configDirNames, configDirNamesGlobal, configFilename,
            inExchangeRate, inStatsRate, mc::SimulationChecking::ExchangeMask(),
            true, restart, restartValue, backupOutput, inMoveStatsVerbosity,
            nbThreads, partitions, partitionsOffset));
  } else if (accepterChoice == "true") {
    return std::unique_ptr<mc::AbstractExchangeSimulation>(
        new MpiBasicExchangeSimulation<mc::OddEvenNeighborLoop, mc::TrueAccept>(
            configDirNames, configDirNamesGlobal, configFilename,
            inExchangeRate, inStatsRate, mc::SimulationChecking::ExchangeMask(),
            true, restart, restartValue, backupOutput, inMoveStatsVerbosity,
            nbThreads, partitions, partitionsOffset));
  } else {
    throw std::invalid_argument(
        fmt::format("{} -- Unknown replex acceptance function. Choose one of "
                    "('hrex', 'remc', 'npt', 'true')",
                    accepterChoice));
  }
}

std::unique_ptr<mc::AbstractExchangeSimulation>
MpiMultidirBuilder(const std::vector<std::string> &configDirNames,
                   const std::vector<std::string> &configDirNamesGlobal,
                   const std::string &configFilename, const bool restart,
                   const std::string &restartValue, const bool backupOutput,
                   const io::MoveStatRecorder::Verbosity inMoveStatsVerbosity,
                   const int nbThreads, const std::vector<int> &partitions,
                   const std::vector<int> &partitionsOffset) {
  // We want to compute all sweep without exchange
  const int exchangeRate = std::numeric_limits<int>::max();
  const int statsRate = std::numeric_limits<int>::max();

  return std::unique_ptr<mc::AbstractExchangeSimulation>(
      new MpiBasicExchangeSimulation<mc::OddEvenNeighborLoop, mc::FalseAccept>(
          configDirNames, configDirNamesGlobal, configFilename, exchangeRate,
          statsRate, mc::SimulationChecking::MultiDirMask(), false, restart,
          restartValue, backupOutput, inMoveStatsVerbosity, nbThreads,
          partitions, partitionsOffset));
}
} // namespace mpi

#endif
