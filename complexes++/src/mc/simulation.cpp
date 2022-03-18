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

#include "mc/simulation.h"
#include "util/timer.h"

namespace mc {

Simulation::Simulation(const std::string &inConfigPath,
                       const std::string &inConfigFilename,
                       const io::MoveStatRecorder::Verbosity inMoveStatsChoice)
    : m_configPath(inConfigPath), m_configFilename(inConfigFilename),
      m_configParameters(
          util::appendPathsIfSubRelative(inConfigPath, m_configFilename)),
      m_shortName(inConfigPath + inConfigFilename),
      m_logFile(util::appendPathsIfSubRelative(
          m_configPath, m_configParameters.value<std::string>("output.log"))),
      m_simuLog(), m_registerLog(RegisterMyLog(m_simuLog)),
      m_sweepTimeCounter(0),
      m_box(io::readBox(util::appendPathsIfSubRelative(
          m_configPath, m_configParameters.value<std::string>("structure")))),
      m_forcefield(io::readForceField(m_configPath, m_configParameters)),
      m_seed(util::fullEntropySeed(
          m_configParameters.value<std::uint32_t>("montecarlo.seed"))),
      isInitialized(false), m_rng(util::RNGEngine(m_seed)),
      m_moveStatsChoice(inMoveStatsChoice) {
  // Print log
  util::Log(util::buildInformation());
  util::Log(util::configInformation(m_configParameters, inConfigFilename));

  util::Log("user seed: {}\n",
            m_configParameters.value<std::uint32_t>("montecarlo.seed"));
  util::Log("program seed: {}\n", m_seed);
}

Simulation::Simulation(io::Deserializer &deserializer)
    : m_configPath(
          deserializer.restore<decltype(m_configPath)>("m_configPath")),
      m_configFilename(
          deserializer.restore<decltype(m_configFilename)>("m_configFilename")),
      m_configParameters(deserializer.restore<decltype(m_configParameters)>(
          "m_configParameters")),
      m_shortName(deserializer.restore<decltype(m_shortName)>("m_shortName")),
      m_logFile(deserializer.restore<decltype(m_logFile)>("m_logFile")),
      m_logBackupName(
          deserializer.restore<decltype(m_logBackupName)>("m_logBackupName")),
      m_registerLog(
          deserializer.restore<decltype(m_registerLog)>("m_registerLog")),
      m_sweepTimeCounter(deserializer.restore<decltype(m_sweepTimeCounter)>(
          "m_sweepTimeCounter")),
      m_box(deserializer.restore<decltype(m_box)>("m_box")),
      m_forcefield(
          deserializer.restore<decltype(m_forcefield)>("m_forcefield")),
      m_seed(deserializer.restore<decltype(m_seed)>("m_seed")),
      m_kernels(deserializer.restore<decltype(m_kernels)>("m_kernels")),
      isInitialized(
          deserializer.restore<decltype(isInitialized)>("isInitialized")),
      m_rng(deserializer.restoreStreamed<decltype(m_rng)>("m_rng")),
      m_moveStatsChoice(deserializer.restore<decltype(m_moveStatsChoice)>(
          "m_moveStatsChoice")) {
  m_mcAlgo = RebuildMcAlgo(deserializer, m_forcefield, m_kernels, m_rng);
}

void Simulation::serialize(io::Serializer &serializer) const {
  serializer.append(m_configPath, "m_configPath");
  serializer.append(m_configFilename, "m_configFilename");
  serializer.append(m_configParameters, "m_configParameters");
  serializer.append(m_shortName, "m_shortName");
  serializer.append(m_logFile, "m_logFile");
  serializer.append(m_logBackupName, "m_logBackupName");
  serializer.append(m_registerLog, "m_registerLog");
  serializer.append(m_sweepTimeCounter, "m_sweepTimeCounter");
  serializer.append(m_box, "m_box");
  serializer.append(m_forcefield, "m_forcefield");
  serializer.append(m_seed, "m_seed");
  serializer.append(m_kernels, "m_kernels");
  serializer.append(isInitialized, "isInitialized");
  serializer.appendStreamed(m_rng, "m_rng");
  serializer.append(m_moveStatsChoice, "m_moveStatsChoice");
  serializer.append(*m_mcAlgo, "m_mcAlgo");
}

Simulation::~Simulation() {
  // Print out move statistics
  if (m_simuLog.isSet()) {
    util::GlobalLog::setGlobalLog(&m_simuLog);

    if (m_mcAlgo != nullptr) {
      util::Log(m_mcAlgo->getMoveStats().getString(m_moveStatsChoice));
    }
  }

  // Print runtime statistics
  if (isInitialized) {
    const auto sweeps_per_min =
        m_mcAlgo->getSessionCounter() / (m_sweepTimeCounter / 60);
    util::Log("Monte Carlo Performance: {} sweep / min\n", sweeps_per_min);
  }
}

int Simulation::rerun() {
  const std::string configFilename = util::appendPathsIfSubRelative(
      m_configPath, m_configParameters.value<std::string>("structure"));
  auto system = io::readTopologies(configFilename, m_forcefield.beadTypes());
  m_kernels = io::readKernelMapping(configFilename, *system.domains());
  mc::McRerun(".", system.domains(), m_box, m_configParameters, m_forcefield,
              m_kernels);
  return 0;
}

int Simulation::run(const int nbSweepToProceed) {
  util::GlobalLog::setGlobalLog(&m_simuLog);

  if (isInitialized == false) {
    throw std::runtime_error(
        "You must init the simulation before calling run.");
  }

  auto timer = util::Timer();
  m_mcAlgo->run(nbSweepToProceed);
  m_sweepTimeCounter += timer.stopAndGetElapsed();

  return 0;
}

//! Init based on a restart file (restartValue)
//! The currentSweepIdx will be shifted accordingly
void Simulation::initFromRestart(const std::string &restartValue) {
  util::GlobalLog::setGlobalLog(&m_simuLog);
  m_simuLog.setLogFile(m_logFile);

  // reset performance counter
  m_sweepTimeCounter = 0;

  if (isInitialized) {
    throw std::runtime_error(
        "Simulation.initFromRestart cannot be called more than once.");
  }

  const std::string configFilename = util::appendPaths(
      m_configPath, m_configParameters.value<std::string>("structure"));

  const auto system =
      io::readTopologies(configFilename, m_forcefield.beadTypes());

  m_kernels = io::readKernelMapping(configFilename, *system.domains());
  // restart file should be in the same directory as the config...
  const std::string restartFilename =
      util::appendPaths(m_configPath, restartValue);
  util::Log("restart-filename: {}\n", restartFilename);

  m_mcAlgo = io::readRestart(restartFilename, m_forcefield, m_kernels, m_rng);

  // we update nstructures to allow extending simulations
  m_mcAlgo->setNStructures(m_configParameters.value<int>("output.nstructures"));
  // Here we must also allow to change our output files. This is useful when one
  // only has the restart file and allows to manually emulate the `-noappend`
  // option of GROMACS.
  const auto traj = m_configParameters.value<std::string>("output.file");
  const auto stat = m_configParameters.value<std::string>("output.stat-file");
  if (traj != m_mcAlgo->getTrajName()) {
    m_mcAlgo->setTraj(traj);
  }
  if (stat != m_mcAlgo->getStatName()) {
    m_mcAlgo->setStat(stat);
  }

  util::Log("restarting simulation at sweep: {}\n",
            m_mcAlgo->getCurrentSweepIdx());
  isInitialized = true;
}

//! Init and backup existing output files if backupOutput is true
void Simulation::init(const bool backupOutput) {
  util::GlobalLog::setGlobalLog(&m_simuLog);

  if (isInitialized) {
    throw std::runtime_error(
        "Simulation.init cannot be called more than once.");
  }

  const std::string configFilename = util::appendPathsIfSubRelative(
      m_configPath, m_configParameters.value<std::string>("structure"));

  auto system = io::readTopologies(configFilename, m_forcefield.beadTypes());

  m_kernels = io::readKernelMapping(configFilename, *system.domains());

  if (backupOutput) {
    util::backUpIfExists(util::appendPathsIfSubRelative(
        m_configPath,
        m_configParameters.value<std::string>("output.stat-file")));
    util::backUpIfExists(util::appendPathsIfSubRelative(
        m_configPath, m_configParameters.value<std::string>("output.file")));
    util::backUpIfExists(m_logFile);
  } else {
    util::Log("overwriting existing files\n");
    util::deleteFileIfExists(util::appendPathsIfSubRelative(
        m_configPath,
        m_configParameters.value<std::string>("output.stat-file")));
    util::deleteFileIfExists(util::appendPathsIfSubRelative(
        m_configPath, m_configParameters.value<std::string>("output.file")));
    if (m_logBackupName.size()) {
      util::deleteFileIfExists(m_logBackupName);
    }
    util::deleteFileIfExists(m_logFile);
  }

  m_simuLog.setLogFile(m_logFile);

  const auto refFileType = m_configParameters.experimental_value<std::string>(
      "output.reference_file_type", "pdb");
  const auto ref_fname =
      util::appendPathsIfSubRelative(
          m_configPath, util::basename(m_configParameters.value<std::string>(
                            "output.file"))) +
      "_reference." + refFileType;
  util::deleteFileIfExists(ref_fname);
  auto refFile = io::TrajectoryFile(ref_fname);
  io::writeModel(refFile, *system.domains(), m_box, 0, 0,
                 m_forcefield.beadTypes());
  util::Log("reference structure saved in: {}\n", ref_fname);
  m_mcAlgo = mc::McBuild(m_configPath, system, m_rng, m_box, m_configParameters,
                         m_forcefield, m_kernels);
  isInitialized = true;
}

const std::string &Simulation::getShortName() const { return m_shortName; }

util::Logger &Simulation::getLogger() { return m_simuLog; }

//! Return the number of domains for this simulation
//! The topology is temporarly loaded if not already in
//! memory
size_t Simulation::getNbDomains() const {
  if (isInitialized == false) {
    // Simulation has not been init
    const auto structure = io::readTopologies(
        util::appendPaths(m_configPath,
                          m_configParameters.value<std::string>("structure")),
        m_forcefield.beadTypes());
    return structure.domains()->size();
  } else {
    return m_mcAlgo->getNbDomains();
  }
}

//! The energy of the system (from the latest run)
//! if no run has been executed DOUBLE_MAX is returned
double Simulation::getEnergy() const { return m_mcAlgo->getEnergy(); }

//! The number of sweeps that the simulation will run
int Simulation::getNbSweep() const { return m_mcAlgo->getNbSweeps(); }

//! Return the number of sweeps in the simulation minus the ones
//! already done if based on a "restart"
int Simulation::getRemainingNbSweep() const {
  return m_mcAlgo->getNbSweeps() - m_mcAlgo->getCurrentSweepIdx();
}

//! Current forcefield
const energy::ForceField &Simulation::getForceField() const {
  return m_mcAlgo->getForceField();
}

//! Compute the energy for the current topology but another
//! force field
energy::rEnergyMatrix
Simulation::computeEnergyForFF(const energy::ForceField &inForceField) const {
  return m_mcAlgo->computeEnergyForFF(inForceField);
}

//! Exchange coordinates (and other needed attributes)
void Simulation::swapCoordinates(Simulation &other,
                                 std::tuple<bool, double> &&isAccepted) {
  std::swap(m_kernels, other.m_kernels);
  m_mcAlgo->swapCoordinates(*other.m_mcAlgo, std::move(isAccepted));
}

//! Exchange coordinates (and other needed attributes)
//! but also provide the energy matrices for each other force field to
//! avoid recomputation
void Simulation::swapCoordinates(
    Simulation &other,
    std::tuple<bool, double, energy::rEnergyMatrix, energy::rEnergyMatrix>
        &&isAcceptedWithEnergy) {
  std::swap(m_kernels, other.m_kernels);
  m_mcAlgo->swapCoordinates(*other.m_mcAlgo, std::move(isAcceptedWithEnergy));
}

//! Get the inverse temperature in reduced units from mc algo (beta)
double Simulation::getBeta() const { return m_mcAlgo->getBeta(); }
util::rvec Simulation::getBox() const { return m_mcAlgo->getBox(); }
double Simulation::getPressure() const { return m_mcAlgo->getPressure(); }

//! Get current random engine
util::RNGEngine &Simulation::getRandEngine() { return m_rng; }
//! Print inStr to current simu log
void Simulation::printToLog(const std::string &inStr) { m_simuLog(inStr); }
} // namespace mc
