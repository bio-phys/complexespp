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
#ifndef MPI_MPIAPPLICATION_H
#define MPI_MPIAPPLICATION_H

#include <fstream>
#include <memory>

// Will include mpi.h
#include "mc/simulation.h"
#include "mc/simulationchecking.h"
#include "mpi/mpicliargs.h"
#include "mpi/mpiexchangebuilder.h"
#include "mpi/mpiutils.h"
#include "parallelization/taskslimiter.h"

namespace mpi {

/**
 * @brief The MpiApplication class is in charge of the high level
 * choice of the simurithm. It decides what kind of simurithm should
 * be executed.
 * @code MpiApplication app(argc, argv);
 * @code return app.run();
 */
class MpiApplication {
  //////////////////////////////////////////////////////////////////////////
  /// Attributes
  //////////////////////////////////////////////////////////////////////////

  const int m_myRank;
  const int m_nbProcesses;

  const MpiCLIArgs m_args;

  //////////////////////////////////////////////////////////////////////////
  /// Private methods
  //////////////////////////////////////////////////////////////////////////

  //! For the remc/exchange executions
  int multiDirExecution() const {
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

    std::vector<int> partitions;
    std::vector<int> partitionsOffset;
    if (m_args.hasKey("mpi-partitions")) {
      partitions = m_args.value<std::vector<int>>("mpi-partitions");
      if (static_cast<int>(partitions.size()) != m_nbProcesses) {
        throw std::runtime_error(
            fmt::format("You are currently executing {} processes and provide "
                        "the partition for {} processes.\n"
                        "You must provide one integer per process (or do not "
                        "use the -partition parameter to use default values).",
                        m_nbProcesses, partitions.size()));
      }

      partitionsOffset.resize(m_nbProcesses + 1);
      partitionsOffset[0] = 0;
      for (int idxProcess = 0; idxProcess < m_nbProcesses; ++idxProcess) {
        partitionsOffset[idxProcess + 1] =
            partitionsOffset[idxProcess] + partitions[idxProcess];
      }

      if (partitionsOffset[m_nbProcesses] != nbSimu) {
        throw std::runtime_error(fmt::format(
            "The number of simulations assigned to the process is {} but it "
            "must be equal to the number of simulations {}.\n"
            "You must adjust your partitions and re-launch the application.",
            partitionsOffset[m_nbProcesses], nbSimu));
      }
    } else {
      partitions.resize(m_nbProcesses);
      partitionsOffset.resize(m_nbProcesses + 1);
      const double partSize =
          static_cast<double>(nbSimu) / static_cast<double>(m_nbProcesses);
      partitionsOffset[0] = 0;
      for (int idxProcess = 0; idxProcess < m_nbProcesses; ++idxProcess) {
        partitionsOffset[idxProcess + 1] =
            static_cast<int>(partSize * static_cast<double>(idxProcess + 1));
        partitions[idxProcess] =
            partitionsOffset[idxProcess + 1] - partitionsOffset[idxProcess];
        DEBUG_ASSERT(partitions[idxProcess] != 0, "Partition cannot be 0");
      }
      DEBUG_ASSERT(partitionsOffset[m_nbProcesses] == nbSimu,
                   "Sum of the partitions must be equal to the number of simu");
    }

    const std::vector<std::string> configDirNamesLocal(
        configDirNames.begin() + partitionsOffset[m_myRank],
        configDirNames.begin() + partitionsOffset[m_myRank + 1]);

    // Build an exchange simu
    std::unique_ptr<mc::AbstractExchangeSimulation> simu;
    if (m_args.hasKey("replex") == false) {
      // pure multidir
      simu = mpi::MpiMultidirBuilder(
          configDirNamesLocal, configDirNames, configFilename,
          m_args.hasKey("restart"),
          m_args.hasKey("restart") ? m_args.value("restart") : "",
          m_args.value<bool>("backup"), moveStatsVerbosity, nbThreads,
          partitions, partitionsOffset);
    } else {
      simu = mpi::MpiExchangeBuilder(
          configDirNamesLocal, configDirNames, configFilename,
          m_args.value<int>("replex"), m_args.value<int>("replex-stat"),
          m_args.value("replex-accept"), m_args.hasKey("restart"),
          m_args.hasKey("restart") ? m_args.value("restart") : "",
          m_args.value<bool>("backup"), moveStatsVerbosity, nbThreads,
          partitions, partitionsOffset);
    }

    // Ensure compatibility
    const bool assertOnFailure = true;
    simu->compareSimulations(assertOnFailure);

    return simu->run();
  }

  //! For the single execution
  int singleSrcExecution() const {
    if (m_nbProcesses != 1) {
      throw std::runtime_error(fmt::format(
          "You are currently executing {} processes for a single simulation.\n"
          "The number of mpi processes cannt be greater than the number of "
          "simulation.",
          m_nbProcesses));
    }

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
      simu.printToLog(
          fmt::format("Parallel -- There are {} threads\n", nbThreads));
      simu.printToLog(fmt::format("Parallel -- There are {} mpi processes\n",
                                  m_nbProcesses));
#pragma omp parallel default(shared) num_threads(nbThreads)
      {
#pragma omp master
        { returnedValue = simu.run(simu.getRemainingNbSweep()); }
      }
      simu.printToLog(util::endingString());
    }
    return returnedValue;
  }

  //////////////////////////////////////////////////////////////////////////
  /// MPI methods to call during the construction
  //////////////////////////////////////////////////////////////////////////

  static int InitMpi(int inArgc, char **inArgv) {
    MPI_ASSERT(MPI_Init(&inArgc, &inArgv));
    int myRank;
    MPI_ASSERT(MPI_Comm_rank(MPI_COMM_WORLD, &myRank));
    return myRank;
  }

  static int GetNbProcesses() {
    int nbProcesses;
    MPI_ASSERT(MPI_Comm_size(MPI_COMM_WORLD, &nbProcesses));
    return nbProcesses;
  }

  //////////////////////////////////////////////////////////////////////////
  /// Public methods
  //////////////////////////////////////////////////////////////////////////
public:
  MpiApplication(int inArgc, char **inArgv, const MPICLIArgs args)
      : m_myRank(InitMpi(inArgc, inArgv)), m_nbProcesses(GetNbProcesses()),
        m_args(args) {}

  ~MpiApplication() {
#ifdef USE_TIMINGOUTPUT
    std::ofstream global_timing_stream("./timings" + std::to_string(m_myRank) +
                                       ".txt");
    GlobalEventManager.show(global_timing_stream);
#endif
    MPI_ASSERT(MPI_Finalize());
  }

  //! Forbid copy
  MpiApplication(const MpiApplication &) = delete;
  MpiApplication &operator=(const MpiApplication &) = delete;

  //! execute the application
  int run() {
    if (m_args.hasKey("version")) {
      if (m_myRank == 0) {
        fmt::print(util::buildInformation());
      }
      return 0;
    }

    if (m_args.hasKey("help") || !m_args.hasKey("config")) {
      if (m_myRank == 0) {
        setup::printHelp(m_args);
      }
      return 0;
    }

    if (m_args.hasKey("multidir")) {
      return multiDirExecution();
    } else {
      return singleSrcExecution();
    }
  }
};
} // namespace mpi

#endif
