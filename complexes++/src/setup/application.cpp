// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------

#include "setup/application.h"
#include "mc/exchangebuilder.h"
#include "mc/simulation.h"
#include "mc/simulationchecking.h"
#include "parallelization/taskslimiter.h"

namespace setup {

int Application::multiDirExecution() const {
  if (m_args.value<bool>("rerun")) {
    throw std::invalid_argument(
        "Rerun option is not supported for multidir.\n");
  }

  const std::string configFilename = m_args.value("config");
  const std::vector<std::string> configDirNames =
      m_args.value<std::vector<std::string>>("multidir");
  const int nbSimu = static_cast<int>(configDirNames.size());
  const int nbThreads = m_args.value<int>("nb-threads");
  util::GlobalLog::setNumberOfThreads(nbThreads);

  if (nbSimu == 0) {
    throw std::invalid_argument(
        "Error in multidir exchange.\n"
        "No directory has been given in parameter.");
  }

  const io::MoveStatRecorder::Verbosity moveStatsVerbosity =
      m_args.getMappingValue<io::MoveStatRecorder::Verbosity>(
          "movestats", {{"none", io::MoveStatRecorder::TXT_NONE},
                        {"perdomain", io::MoveStatRecorder::TXT_PER_DOMAIN},
                        {"pertype", io::MoveStatRecorder::TXT_PER_TYPE},
                        {"all", io::MoveStatRecorder::TXT_ALL}});

  // Build an exchange algo
  std::unique_ptr<mc::AbstractExchangeSimulation> simu;
  if (m_args.hasKey("replex") == false) {
    // pure multidir
    simu = mc::MultidirBuilder(
        configDirNames, configFilename, m_args.hasKey("restart"),
        m_args.hasKey("restart") ? m_args.value("restart") : "",
        m_args.value<bool>("backup"), moveStatsVerbosity, nbThreads);
  } else {
    simu = mc::ExchangeBuilder(
        configDirNames, configFilename, m_args.value<int>("replex"),
        m_args.value<int>("replex-stat"), m_args.value("replex-accept"),
        m_args.hasKey("restart"),
        m_args.hasKey("restart") ? m_args.value("restart") : "",
        m_args.value<bool>("backup"), moveStatsVerbosity, nbThreads,
        m_args.value<std::string>("replex-verbosity"));
  }

  // Ensure compatibility
  const bool assertOnFailure = true;
  simu->compareSimulations(assertOnFailure);

  return simu->run();
}

int Application::singleSrcExecution() const {
  const io::MoveStatRecorder::Verbosity moveStatsVerbosity =
      m_args.getMappingValue<io::MoveStatRecorder::Verbosity>(
          "movestats", {{"none", io::MoveStatRecorder::TXT_NONE},
                        {"perdomain", io::MoveStatRecorder::TXT_PER_DOMAIN},
                        {"pertype", io::MoveStatRecorder::TXT_PER_TYPE},
                        {"all", io::MoveStatRecorder::TXT_ALL}});
  mc::Simulation simu(util::dirname(m_args.value("config")),
                      util::filename(m_args.value("config")),
                      moveStatsVerbosity);

  int returnedValue;
  const int nbThreads = m_args.value<int>("nb-threads");
  util::GlobalLog::setNumberOfThreads(nbThreads);

  // Init TIMEZONE before a parallel section
  TIMEZONE_OMP_INIT_PREPARALLEL(nbThreads)
  // Force task creation in parallel section if more than one thread
  vectorization::TasksLimiter::Controller.setEnableTasks(nbThreads > 1);

  if (m_args.value<bool>("rerun")) {
#pragma omp parallel default(shared) num_threads(nbThreads)
    {
#pragma omp master
      { returnedValue = simu.rerun(); }
    }
  } else {
    if (m_args.hasKey("restart")) {
      simu.initFromRestart(m_args.value("restart"));
    } else {
      simu.init(m_args.value<bool>("backup"));
    }
    simu.printToLog(util::startingString());
#pragma omp parallel default(shared) num_threads(nbThreads)
    {
#pragma omp master
      { returnedValue = simu.run(simu.getRemainingNbSweep()); }
    }
    simu.printToLog(util::endingString());
  }
  return returnedValue;
}

int Application::run() {
  if (m_args.hasKey("version")) {
    fmt::print(util::buildInformation());
    return 0;
  }

  if (m_args.hasKey("help") || !m_args.hasKey("config")) {
    printHelp(m_args);
    return 0;
  }

  if (m_args.hasKey("multidir")) {
    return multiDirExecution();
  } else {
    return singleSrcExecution();
  }
}
}  // namespace setup
