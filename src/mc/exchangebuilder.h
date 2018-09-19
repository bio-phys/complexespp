// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#ifndef EXCHANGEBUILDER_H
#define EXCHANGEBUILDER_H

#include <memory>
#include <string>
#include <vector>

#include "mc/abstractexchangesimulation.h"
#include "mc/basicexchangesimulation.h"
#include "mc/exchangeaccepter.h"
#include "mc/exchangelooper.h"
#include "mc/simulationchecking.h"
#include "setup/cliargs.h"

namespace mc {

std::unique_ptr<AbstractExchangeSimulation> ExchangeBuilder(
    const std::vector<std::string>& configDirNames,
    const std::string& configFilename, const int inExchangeRate,
    const int inStatsRate, const std::string& accepterChoice,
    const bool restart, const std::string& restartValue,
    const bool backupOutput,
    const io::MoveStatRecorder::Verbosity inMoveStatsVerbosity,
    const int nbThreads, const std::string& logVerbosity) {
  if (accepterChoice == "hrex") {
    return std::unique_ptr<AbstractExchangeSimulation>(
        new BasicExchangeSimulation<OddEvenNeighborLoop, HREXAccept>(
            configDirNames, configFilename, inExchangeRate, inStatsRate,
            SimulationChecking::ExchangeMask(), logVerbosity, restart,
            restartValue, backupOutput, inMoveStatsVerbosity, nbThreads));
  } else if (accepterChoice == "remc") {
    return std::unique_ptr<AbstractExchangeSimulation>(
        new BasicExchangeSimulation<OddEvenNeighborLoop, REMCAccept>(
            configDirNames, configFilename, inExchangeRate, inStatsRate,
            SimulationChecking::ExchangeMask(), logVerbosity, restart,
            restartValue, backupOutput, inMoveStatsVerbosity, nbThreads));
  } else if (accepterChoice == "npt") {
    return std::unique_ptr<AbstractExchangeSimulation>(
        new BasicExchangeSimulation<OddEvenNeighborLoop, NPTAccept>(
            configDirNames, configFilename, inExchangeRate, inStatsRate,
            SimulationChecking::ExchangeMask(), logVerbosity, restart,
            restartValue, backupOutput, inMoveStatsVerbosity, nbThreads));
  } else if (accepterChoice == "true") {
    return std::unique_ptr<AbstractExchangeSimulation>(
        new BasicExchangeSimulation<OddEvenNeighborLoop, TrueAccept>(
            configDirNames, configFilename, inExchangeRate, inStatsRate,
            SimulationChecking::ExchangeMask(), logVerbosity, restart,
            restartValue, backupOutput, inMoveStatsVerbosity, nbThreads));
  } else {
    throw std::invalid_argument(
        fmt::format("{} -- Unknown replex acceptance function. Choose one of "
                    "('hrex', 'remc', 'npt', 'true')",
                    accepterChoice));
  }
}

std::unique_ptr<AbstractExchangeSimulation> MultidirBuilder(
    const std::vector<std::string>& configDirNames,
    const std::string& configFilename, const bool restart,
    const std::string& restartValue, const bool backupOutput,
    const io::MoveStatRecorder::Verbosity inMoveStatsVerbosity,
    const int nbThreads) {
  // We want to compute all sweep without exchange
  const int exchangeRate = std::numeric_limits<int>::max();
  const int statsRate = std::numeric_limits<int>::max();

  return std::unique_ptr<AbstractExchangeSimulation>(
      new BasicExchangeSimulation<OddEvenNeighborLoop, FalseAccept>(
          configDirNames, configFilename, exchangeRate, statsRate,
          SimulationChecking::MultiDirMask(), "none", restart, restartValue,
          backupOutput, inMoveStatsVerbosity, nbThreads));
}
}  // namespace mc

#endif
