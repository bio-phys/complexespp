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
#ifndef BASICEXCHANGEALGORITHM_H
#define BASICEXCHANGEALGORITHM_H

#include <string>
#include <vector>

#include "io/movestatrecorder.h"
#include "mc/abstractexchangesimulation.h"
#include "mc/exchangelogger.h"
#include "mc/simulation.h"
#include "mc/simulationchecking.h"
#include "parallelization/taskslimiter.h"

namespace mc {

enum class ExchangeLogVerbosity { all, stats, none };

ExchangeLogVerbosity determineVerbosity(const std::string &val) {
  if (val == "stats") {
    return ExchangeLogVerbosity::stats;
  } else if (val == "all") {
    return ExchangeLogVerbosity::all;
  } else if (val == "none") {
    return ExchangeLogVerbosity::none;
  } else {
    throw std::invalid_argument(
        "Unknown verbosty level. Choose one of 'stats, all, none'");
  }
}

//! This class lets compute an exchanged based algorithm
//! with templatized exchange iterations and accept function.
template <class LoopClass, class AcceptClass>
class BasicExchangeSimulation : public AbstractExchangeSimulation {
  //! The number of simulations
  const int m_nbSimu;
  //! The simulations
  std::vector<Simulation> m_simus;
  //! The exchange rate (number of sweeps between exchange)
  const int m_exchangeRate;
  //! The stats rate (number of sweeps between stats output)
  const int m_statisticRate;
  //! The simulation directories
  const std::vector<std::string> m_configDirNames;
  //! The config filename
  const std::string m_configFilename;
  //! The accept object
  AcceptClass m_accepter;
  //! The given comparison mask
  const SimulationChecking::ComparisonFlag m_defaultComparisonMask;
  //! To know if we should ouput exchange log
  const ExchangeLogVerbosity m_enableExchangeLog;
  //! The number of threads we can create
  const int m_nbThreads;

public:
  BasicExchangeSimulation(
      const std::vector<std::string> &configDirNames,
      const std::string &configFilename, const int inExchangeRate,
      const int inStatsRate,
      const SimulationChecking::ComparisonFlag inComparisonMask,
      const std::string inEnableExchangeLog, const bool restart,
      const std::string &restartValue, const bool backupOutput,
      const io::MoveStatRecorder::Verbosity inMoveStatsVerbosity,
      const int inNbThreads)
      : m_nbSimu(static_cast<int>(configDirNames.size())),
        m_exchangeRate(inExchangeRate), m_statisticRate(inStatsRate),
        m_configDirNames(configDirNames), m_configFilename(configFilename),
        m_defaultComparisonMask(inComparisonMask),
        m_enableExchangeLog(determineVerbosity(inEnableExchangeLog)),
        m_nbThreads(inNbThreads) {
    TIMEZONE("BasicExchangeSimulation")
    // Build the simulation
    m_simus.reserve(m_nbSimu);

    for (int idxSimu = 0; idxSimu < m_nbSimu; ++idxSimu) {
      m_simus.emplace_back(configDirNames[idxSimu], configFilename,
                           inMoveStatsVerbosity);
    }

    for (int idxSimu = 0; idxSimu < m_nbSimu; idxSimu += 1) {
      TIMEZONE(fmt::format("simu_init {}", m_simus[idxSimu].getShortName()))
      util::GlobalLog::setGlobalLog(&m_simus[idxSimu].getLogger());
      if (restart) {
        m_simus[idxSimu].initFromRestart(restartValue);
      } else {
        m_simus[idxSimu].init(backupOutput);
      }
    }
  }

  ~BasicExchangeSimulation() {}

  bool compareSimulations(const bool assertOnFailure) final {
    // we delegate the comparison to the CompareExchange function
    auto comparisonMask = m_defaultComparisonMask;
    if (m_accepter.needSameForceField()) {
      comparisonMask =
          comparisonMask | SimulationChecking::ComparisonFlag::FORCEFIELD;
    }
    return SimulationChecking::compareSimulations(m_simus, comparisonMask,
                                                  assertOnFailure);
  }

  int run() final {
    TIMEZONE("run")
    //! The output log system
    ExchangeLogger logger(m_nbSimu);
    if (m_enableExchangeLog != ExchangeLogVerbosity::none) {
      for (int idxSimu = 0; idxSimu < m_nbSimu; ++idxSimu) {
        // Tell the log system where to output the current simu output
        logger.setLoggerRef(idxSimu, m_simus[idxSimu].getLogger());
      }
    }

    for (int idxSimu = 0; idxSimu < m_nbSimu; ++idxSimu) {
      m_simus[idxSimu].printToLog(util::startingString());
    }

    // Each simulations should have the same number of sweeps
    const int nbSweeps = m_simus[0].getNbSweep();

    if (m_enableExchangeLog != ExchangeLogVerbosity::none) {
      // Print out the log header
      logger.printHeader(m_configDirNames, m_configFilename, m_exchangeRate,
                         m_statisticRate, nbSweeps);
    }

    // Compute nbSweeps sweeps by step m_exchangeRate
    for (int idxSweep = 0; idxSweep < nbSweeps; idxSweep += m_exchangeRate) {
      const int nbSweepToProceed =
          std::min(m_exchangeRate, nbSweeps - idxSweep);

      // Init TIMEZONE before a parallel section
      TIMEZONE_OMP_INIT_PREPARALLEL(m_nbThreads)
      // Force task creation in parallel section if more than one thread per
      // simu
      vectorization::TasksLimiter::Controller.setEnableTasks(m_nbThreads >
                                                             m_nbSimu);
// Compute the mc sweep
#pragma omp parallel num_threads(m_nbThreads) default(shared)
      {
#pragma omp for schedule(dynamic, 1) nowait
        for (int idxSimu = 0; idxSimu < m_nbSimu; idxSimu += 1) {
          TIMEZONE(fmt::format("run {}", m_simus[idxSimu].getShortName()))
          util::GlobalLog::setGlobalLog(&m_simus[idxSimu].getLogger());
          m_simus[idxSimu].run(nbSweepToProceed);
        }
        // When a thread becomes available we turn the tasks to ON
        vectorization::TasksLimiter::Controller.setEnableTasks(true);
      }

      // Exchange only if not the last loop
      if (idxSweep + nbSweepToProceed != nbSweeps) {
        if (m_enableExchangeLog != ExchangeLogVerbosity::none) {
          logger.startExchange(idxSweep + nbSweepToProceed);
        }
        TIMEZONE("exchange")
        // Iterate over the simulations
        for (auto iterExchange :
             LoopClass(m_nbSimu, idxSweep / m_exchangeRate)) {
          const int idx1 = iterExchange.first;
          const int idx2 = iterExchange.second;
          // Ask if we need to exchange
          auto resExchange = m_accepter.accept(m_simus[idx1], m_simus[idx2],
                                               m_simus[idx1].getRandEngine());
          const bool doExchange = std::get<0>(resExchange);
          const double probability = std::get<1>(resExchange);
          if (m_enableExchangeLog != ExchangeLogVerbosity::none) {
            // Print the result in the log
            logger.addAttempt(idx1, idx2, m_simus[idx1].getEnergy(),
                              m_simus[idx2].getEnergy(), probability,
                              doExchange);
          }
          // Perform the exchange if needed
          if (doExchange) {
            m_simus[idx1].swapCoordinates(m_simus[idx2],
                                          std::move(resExchange));
          }
        }

        if (m_enableExchangeLog != ExchangeLogVerbosity::none) {
          // Print the attempts
          if (m_enableExchangeLog == ExchangeLogVerbosity::all) {
            logger.printAttempts(LoopClass::IsOddEvenLoop());
          }
          if (m_enableExchangeLog == ExchangeLogVerbosity::all ||
              m_enableExchangeLog == ExchangeLogVerbosity::stats)
            if (((idxSweep + 1) % m_statisticRate) == 0) {
              // Print the stats
              logger.printStats();
            }
        }
      }
    }

    if (m_enableExchangeLog != ExchangeLogVerbosity::none) {
      logger.printStats();
    }

    for (int idxSimu = 0; idxSimu < m_nbSimu; ++idxSimu) {
      m_simus[idxSimu].printToLog(util::endingString());
    }

    return 0;
  }
};
} // namespace mc

#endif
