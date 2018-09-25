// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#ifndef MPIBASICEXCHANGESIMULATION_H
#define MPIBASICEXCHANGESIMULATION_H

#include <string>
#include <vector>

#include "io/serializer.h"
#include "mc/abstractexchangesimulation.h"
#include "mc/exchangelogger.h"
#include "mc/simulation.h"
#include "mc/simulationchecking.h"
#include "mpi/mpiutils.h"
#include "parallelization/taskslimiter.h"

namespace mpi {

//! This class lets compute an exchanged based algorithm
//! with templatized exchange iterations and accept function.
template <class LoopClass, class AcceptClass>
class MpiBasicExchangeSimulation : public mc::AbstractExchangeSimulation {
  //! The number of simulations
  const int m_nbSimu;
  //! The simulations
  std::vector<mc::Simulation> m_simus;
  //! The exchange rate (number of sweeps between exchange)
  const int m_exchangeRate;
  //! The stats rate (number of sweeps between stats output)
  const int m_statisticRate;
  //! The simulation directories
  const std::vector<std::string> m_configDirNames;
  //! The simulation directories for all ranks
  const std::vector<std::string> m_configDirNamesGlobal;
  //! The config filename
  const std::string m_configFilename;
  //! The accept object
  AcceptClass m_accepter;
  //! The given comparison mask
  const mc::SimulationChecking::ComparisonFlag m_defaultComparisonMask;
  //! To know if we should ouput exchange log
  const bool m_enableExchangeLog;
  //! The number of threads we can create
  const int m_nbThreads;

  //! Number of simu per process
  const std::vector<int> m_partitions;
  //! Offset of number of simu per process
  //! (m_partitions[i] = m_partitionsOffset[i+1]-m_partitionsOffset[i])
  const std::vector<int> m_partitionsOffset;
  //! Total number of simulations for all processes
  const int m_nbSimuGlobal;
  //! My mpi rank
  int m_myRank;
  //! Nb mpi processes
  int m_nbProcesses;

  void coreExchange(mc::ExchangeLogger& loger, const int idx1, const int idx2,
                    mc::Simulation& simu1, mc::Simulation& simu2) {
    // Ask if we need to exchange
    auto resExchange = m_accepter.accept(simu1, simu2, simu1.getRandEngine());
    const bool doExchange = std::get<0>(resExchange);
    const double probability = std::get<1>(resExchange);
    if (m_enableExchangeLog) {
      // Print the result in the log
      loger.addAttempt(idx1, idx2, simu1.getEnergy(), simu2.getEnergy(),
                       probability, doExchange);
    }
    // Perform the exchange if needed
    if (doExchange) {
      simu1.swapCoordinates(simu2, std::move(resExchange));
    }
  }

  mc::Simulation receiveSimulation(const int idxProcessSrc) const {
    TIMEZONE("receiveSimulation")
    MPI_Status status;
    MPI_ASSERT(MPI_Probe(idxProcessSrc, 1, MPI_COMM_WORLD, &status));

    int sizeToRecv;
    MPI_ASSERT(MPI_Get_count(&status, MPI_UNSIGNED_CHAR, &sizeToRecv));

    std::vector<unsigned char> buffer;
    buffer.resize(sizeToRecv);

    {
      TIMEZONE("MPI_Recv")
      MPI_ASSERT(MPI_Recv(buffer.data(), sizeToRecv, MPI_UNSIGNED_CHAR,
                          idxProcessSrc, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE));
    }

    TIMEZONE("deserializer")
    io::Deserializer deserializer(buffer.data(), buffer.size());
    return deserializer.restore<mc::Simulation>("simulation");
  }

  void sendSimulation(const mc::Simulation& simuToSend,
                      const int idxProcessDest) const {
    TIMEZONE("sendSimulation")
    io::Serializer serializer;
    serializer.append(simuToSend, "simulation");
    const std::vector<unsigned char> buffer = serializer.releaseBuffer();

    DEBUG_ASSERT(buffer.size() < std::numeric_limits<int>::max(),
                 "Invalid size of buffer");

    TIMEZONE("MPI_Send")
    MPI_ASSERT(MPI_Send(const_cast<unsigned char*>(buffer.data()),
                        static_cast<int>(buffer.size()), MPI_UNSIGNED_CHAR,
                        idxProcessDest, 1, MPI_COMM_WORLD));
  }

 public:
  MpiBasicExchangeSimulation(
      const std::vector<std::string>& configDirNames,
      const std::vector<std::string>& configDirNamesGlobal,
      const std::string& configFilename, const int inExchangeRate,
      const int inStatsRate,
      const mc::SimulationChecking::ComparisonFlag inComparisonMask,
      const bool inEnableExchangeLog, const bool restart,
      const std::string& restartValue, const bool backupOutput,
      const io::MoveStatRecorder::Verbosity inMoveStatsVerbosity,
      const int inNbThreads, const std::vector<int>& inPartitions,
      const std::vector<int>& inPartitionsOffset)
      : m_nbSimu(static_cast<int>(configDirNames.size())),
        m_exchangeRate(inExchangeRate),
        m_statisticRate(inStatsRate),
        m_configDirNames(configDirNames),
        m_configDirNamesGlobal(configDirNamesGlobal),
        m_configFilename(configFilename),
        m_defaultComparisonMask(inComparisonMask),
        m_enableExchangeLog(inEnableExchangeLog),
        m_nbThreads(inNbThreads),
        m_partitions(inPartitions),
        m_partitionsOffset(inPartitionsOffset),
        m_nbSimuGlobal(m_partitionsOffset.back()),
        m_myRank(-1) {
    TIMEZONE("MpiBasicExchangeSimulation");
    MPI_ASSERT(MPI_Comm_rank(MPI_COMM_WORLD, &m_myRank));
    MPI_ASSERT(MPI_Comm_size(MPI_COMM_WORLD, &m_nbProcesses));
    // Build the simulation
    m_simus.reserve(m_nbSimu);

    for (int idxSimu = 0; idxSimu < m_nbSimu; ++idxSimu) {
      m_simus.emplace_back(configDirNames[idxSimu], configFilename,
                           inMoveStatsVerbosity);
    }

    for (int idxSimu = 0; idxSimu < m_nbSimu; idxSimu += 1) {
      TIMEZONE(fmt::format("simu_init {}", m_simus[idxSimu].getShortName()));
      util::GlobalLog::setGlobalLog(&m_simus[idxSimu].getLogger());
      if (restart) {
        m_simus[idxSimu].initFromRestart(restartValue);
      } else {
        m_simus[idxSimu].init(backupOutput);
      }
    }
  }

  ~MpiBasicExchangeSimulation() {}

  bool compareSimulations(const bool assertOnFailure) final {
    // we delegate the comparison to the CompareExchange function
    mc::SimulationChecking::ComparisonFlag comparisonMask =
        m_defaultComparisonMask;
    if (m_accepter.needSameForceField()) {
      comparisonMask =
          comparisonMask | mc::SimulationChecking::ComparisonFlag::FORCEFIELD;
    }
    return mc::SimulationChecking::compareSimulations(m_simus, comparisonMask,
                                                      assertOnFailure);
  }

  //! Return true if the simulation of index inIdxSimu
  //! is assigned to this object
  bool simulationIsMine(const int inIdxSimu) const {
    return m_partitionsOffset[m_myRank] <= inIdxSimu &&
           inIdxSimu < m_partitionsOffset[m_myRank + 1];
  }

  //! Return the process id owner of the simu of index inIdxSimu
  int simulationOwner(const int inIdxSimu) const {
    for (int idxProc = 0; idxProc < static_cast<int>(m_partitions.size());
         ++idxProc) {
      if (m_partitionsOffset[idxProc] <= inIdxSimu &&
          inIdxSimu < m_partitionsOffset[idxProc + 1]) {
        return idxProc;
      }
    }
    return -1;
  }

  int run() final {
    TIMEZONE("run")
    //! The output log system
    mc::ExchangeLogger loger(m_nbSimuGlobal);
    if (m_enableExchangeLog) {
      for (int idxSimu = 0; idxSimu < m_nbSimu; ++idxSimu) {
        // Tell the log system where to output the current simu output
        loger.setLoggerRef(idxSimu + m_partitionsOffset[m_myRank],
                           m_simus[idxSimu].getLogger());
      }
    }

    for (int idxSimu = 0; idxSimu < m_nbSimu; ++idxSimu) {
      m_simus[idxSimu].printToLog(util::startingString());
      m_simus[idxSimu].printToLog(
          fmt::format("Parallel -- There are {} threads\n", m_nbThreads));
      m_simus[idxSimu].printToLog(fmt::format(
          "Parallel -- There are {} mpi processes\n", m_nbProcesses));
    }

    // Each simulations should have the same number of sweeps
    const int nbSweeps = m_simus[0].getNbSweep();

    if (m_enableExchangeLog) {
      // Print out the log header
      loger.printHeader(m_configDirNamesGlobal, m_configFilename,
                        m_exchangeRate, m_statisticRate, nbSweeps);
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
        if (m_enableExchangeLog) {
          loger.startExchange(idxSweep);
        }
        TIMEZONE("exchange")
        // Iterate over the simulations
        for (auto iterExchange :
             LoopClass(m_nbSimuGlobal, idxSweep / m_exchangeRate)) {
          const int idx1 = iterExchange.first;
          const int idx2 = iterExchange.second;
          if (simulationIsMine(idx1) && simulationIsMine(idx2)) {
            // Both simulation are mine, proceed as usual
            coreExchange(
                loger, idx1, idx2,
                m_simus[idx1 - m_partitionsOffset[m_myRank]],   // by reference
                                                                // but will be
                                                                // modified
                m_simus[idx2 - m_partitionsOffset[m_myRank]]);  // by reference
                                                                // but will be
                                                                // modified
          } else if (simulationIsMine(idx1)) {
            sendSimulation(m_simus[idx1 - m_partitionsOffset[m_myRank]],
                           simulationOwner(idx2));

            mc::Simulation receivedSimulation =
                receiveSimulation(simulationOwner(idx2));
            // Only one of them is mine
            coreExchange(
                loger, idx1, idx2,
                m_simus[idx1 - m_partitionsOffset[m_myRank]],  // by reference
                                                               // but will be
                                                               // modified
                receivedSimulation);  // by reference but will be modified
          } else if (simulationIsMine(idx2)) {
            mc::Simulation receivedSimulation =
                receiveSimulation(simulationOwner(idx1));
            sendSimulation(m_simus[idx2 - m_partitionsOffset[m_myRank]],
                           simulationOwner(idx1));
            // Only one of them is mine
            coreExchange(
                loger, idx1, idx2,
                receivedSimulation,  // by reference but will be modified
                m_simus[idx2 - m_partitionsOffset[m_myRank]]);  // by reference
                                                                // but will be
                                                                // modified
          }
        }

        if (m_enableExchangeLog) {
          // Print the attempts
          loger.printAttempts(LoopClass::IsOddEvenLoop());
          if (((idxSweep + 1) % m_statisticRate) == 0) {
            // Print the stats
            loger.printStats();
          }
        }
      }
    }

    if (m_enableExchangeLog) {
      loger.printStats();
    }

    for (int idxSimu = 0; idxSimu < m_nbSimu; ++idxSimu) {
      m_simus[idxSimu].printToLog(util::endingString());
    }

    return 0;
  }
};
}  // namespace mpi

#endif  // MPI_BASICEXCHANGEALGORITHM_h
