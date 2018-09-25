// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#ifndef MC_SIMULATION_H
#define MC_SIMULATION_H

#include "domains/abstractdomain.h"
#include "io/io.h"
#include "io/serializer.h"
#include "mc/mc.h"
#include "setup/cliargs.h"
#include "setup/config.h"
#include "util/file.h"
#include "util/log.h"
#include "util/random.h"
#include "version.h"

// For the build TODO !
#include "mc/accept_func.h"
#include "mc/nvt.h"

namespace mc {

/**
 * @brief The Simulation class executes a simulation
 * from the high level view (that is for one config file
 * and its related structure files).
 *
 * The topology loaded during the init stage (not in the
 * constructor).
 * Therefore a choice must be done to select a init
 * from restart or not.
 * Only one call to init is possible, a second one will
 * throw an exception.
 *
 * After a init, it is possible to call "run" several times.
 *
 * It is the responsability to the caller to be in the right
 * directory when creating a Simulation object.
 * However, all other methods of Simulation perform the
 * directory change if needed (with a rollback).
 */
class Simulation : public io::AbstractSerializable {
 protected:
  //! Where the config files are located
  const std::string m_configPath;
  //! The name of the config file
  const std::string m_configFilename;
  //! The config file loaded
  const setup::Config m_configParameters;
  //! A string to refer to the current simu
  //! (it is unique for now but it is not guarantee to remain)
  const std::string m_shortName;

  //! The output log file (given from config)
  const std::string m_logFile;
  //! The name of the backedup log file
  const std::string m_logBackupName;
  //! Log class for this simulation
  util::Logger m_simuLog;
  //! To register log in the constructor
  bool m_registerLog;
  //! To accumulate runtime of sweeps
  double m_sweepTimeCounter;

  //! The size of the simulation box (given from input or restart-file)
  util::rvec m_box;
  //! The forcefield (given from config)
  const energy::ForceField m_forcefield;
  //! The random seed (given from config)
  const std::uint32_t m_seed;
  //! The kernel pairs manager
  pairkernels::PairKernelManager m_kernels;

  //! The underlying MC algo
  std::unique_ptr<mc::AbstractMcAlgo> m_mcAlgo;
  //! To know if the current class has been init
  bool isInitialized;
  //! Random generator initialized from m_seed
  util::RNGEngine m_rng;

  //! The move stats to print
  const io::MoveStatRecorder::Verbosity m_moveStatsChoice;

  static bool RegisterMyLog(util::Logger& currentLog) {
    util::GlobalLog::setGlobalLog(&currentLog);
    return true;
  }

 public:
  //! Pre-load the config file inConfigFilename
  //! If inChangeLog is true (default) then the log is redirected to
  //! config.output.log
  Simulation(const std::string& inConfigPath,
             const std::string& inConfigFilename,
             const io::MoveStatRecorder::Verbosity inMoveStatsChoice);
  ~Simulation();

  //! Forbid copy
  Simulation(const Simulation&) = delete;
  Simulation& operator=(const Simulation&) = delete;
  //! Agree to move
  Simulation(Simulation&&) = default;
  Simulation& operator=(Simulation&&) = default;

  /////////////////////////////////////////////////////////////////
  /// Getter/setter
  /////////////////////////////////////////////////////////////////

  const std::string& getShortName() const;

  util::Logger& getLogger();

  //! Return the number of domains for this simulation
  //! The topology is temporarly loaded if not already in
  //! memory
  size_t getNbDomains() const;

  //! The energy of the system (from the latest run)
  //! if no run has been executed DOUBLE_MAX is returned
  double getEnergy() const;

  //! The number of sweeps that the simulation will run
  int getNbSweep() const;

  //! Return the number of sweeps in the simulation minus the ones
  //! already done if based on a "restart"
  int getRemainingNbSweep() const;

  //! Current forcefield
  const energy::ForceField& getForceField() const;

  //! Compute the energy for the current topology but another
  //! force field
  energy::rEnergyMatrix computeEnergyForFF(
      const energy::ForceField& inForceField) const;

  //! Exchange coordinates (and other needed attributes)
  void swapCoordinates(Simulation& other,
                       std::tuple<bool, double>&& isAccepted);

  //! Exchange coordinates (and other needed attributes)
  //! but also provide the energy matrices for each other force field to
  //! avoid recomputation
  void swapCoordinates(
      Simulation& other,
      std::tuple<bool, double, energy::rEnergyMatrix, energy::rEnergyMatrix>&&
          isAcceptedWithEnergy);

  //! Get the inverse temperature in reduced units from mc algo (beta)
  double getBeta() const;
  util::rvec getBox() const;
  double getPressure() const;

  //! Get current random engine
  util::RNGEngine& getRandEngine();
  //! Print inStr to current simu log
  void printToLog(const std::string& inStr);

  /////////////////////////////////////////////////////////////////
  /// Init
  /////////////////////////////////////////////////////////////////

  //! Init based on a restart file (restartValue)
  //! The currentSweepIdx will be shifted accordingly
  void initFromRestart(const std::string& restartValue);
  void init(const bool backupOutput);

  /////////////////////////////////////////////////////////////////
  /// Run
  /////////////////////////////////////////////////////////////////

  int run(const int nbSweepToProceed);
  int rerun();

  /////////////////////////////////////////////////////////////////
  /// Serialize
  /////////////////////////////////////////////////////////////////
  void serialize(io::Serializer& serializer) const final;
  Simulation(io::Deserializer& deserializer);

  static std::unique_ptr<mc::AbstractMcAlgo> RebuildMcAlgo(
      io::Deserializer& deserializer, const energy::ForceField& inForcefield,
      const pairkernels::PairKernelManager& inKernels, util::RNGEngine& inRng) {
    deserializer.access("m_mcAlgo");
    return mc::AbstractMcAlgo::Rebuild(deserializer, inForcefield, inKernels,
                                       inRng);
  }
};

}  // namespace mc

#endif  // MC_SIMULATION
