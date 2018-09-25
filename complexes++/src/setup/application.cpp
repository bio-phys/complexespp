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

#include "setup/application.h"
#include "mc/exchangebuilder.h"
#include "mc/simulation.h"
#include "mc/simulationchecking.h"
#include "parallelization/taskslimiter.h"

namespace setup {

void printGPL() {
  std::cout << " Copyright (c) 2018 the complexes++ development team and "
               "contributors \n";
  std::cout << " (see the file AUTHORS for the full list of names) \n";
  std::cout << " \n";
  std::cout << " complexes++ is free software: you can redistribute it and/or "
               "modify \n";
  std::cout << " it under the terms of the Lesser GNU General Public License "
               "as published by \n";
  std::cout << " the Free Software Foundation, either version 3 of the "
               "License, or \n";
  std::cout << " (at your option) any later version. \n";
  std::cout << " \n";
  std::cout
      << " complexes++ is distributed in the hope that it will be useful, \n";
  std::cout
      << " but WITHOUT ANY WARRANTY; without even the implied warranty of \n";
  std::cout
      << " MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the \n";
  std::cout << " GNU General Public License for more details. \n";
  std::cout << " \n";
  std::cout << " You should have received a copy of the GNU General Public "
               "License \n";
  std::cout << " along with complexes++.  If not, see "
               "<https://www.gnu.org/licenses/> \n";
  std::cout << std::endl;
}

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
    throw std::invalid_argument("Error in multidir exchange.\n"
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
  printGPL();

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
} // namespace setup
