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
#ifndef ABSTRACTMCALGO_H
#define ABSTRACTMCALGO_H

#include "constants.h"
#include "domains/abstractdomain.h"
#include "energy/energymatrix.h"
#include "energy/forcefield.h"
#include "energy/pairparameter.h"
#include "interactions/interactionbuilder.h"
#include "io/io.h"
#include "io/movestatrecorder.h"
#include "io/serializer.h"
#include "pairkernels/pairkernelmanager.h"
#include "setup/config.h"
#include "util/random.h"

// For the build function
#include "domains/rigiddomain.h"
#include "io/rebuilder.h"

namespace mc {

//! The abstractMcAlgo class is the core part of the Mc algorithm
//! it has a sweep method where the inner exchange are done.
//! To execute it one must call the run method and ask for a given
//! number of sweep.
//! However, this class should mostly be executed by the Simulation class.
//!
//! The accept function follows the template method pattern
//! and must be provided by a subclass.
//!
//! AbstractMcAlgo proposes a swapCoordinates which should be used
//! by exchange algorithms to swap the coordinates.
//! The first call to run (and in case updates have changed the
//! attributes like after an swapCoordinates) a computeAll is executed.
//!
//! In the constructor, the parameter inStartingSweepIdxRef indicates
//! if we start a new simulation or restart.
//! If it is equal to 0, some more inits are done, and the current energy
//! is computed in the constructor to be written as energy for sweep 0.
//! The index is then inc to 1 as first real sweep computation.
//! If the value of inStartingSweepIdxRef is >= 1 then it means
//! we restart a simulation and so no headers are written in the stat file
//! because it should have been already the case.
class AbstractMcAlgo
    : public io::RebuilderCore<AbstractMcAlgo, const energy::ForceField &,
                               const pairkernels::PairKernelManager &,
                               util::RNGEngine &>,
      public io::AbstractSerializable {
protected:
  //! The current topology
  std::shared_ptr<domains::Domains> m_doms;
  std::vector<domains::Topology> m_topologies;
  //! The current energy between the domains of the topology
  //! for the current forcefield
  energy::rEnergyMatrix m_energy;
  //! The simulation box
  util::rvec m_box;
  //! The current force field
  const energy::ForceField &m_forcefield;
  //! The random engine
  util::RNGEngine &m_rng;
  //! The interaction computer for the given topology
  std::unique_ptr<AbstractInteractionAlgorithm<double>> m_interactionComputer;

  //! To know if a computeAll must be done in run
  bool m_needComputeAll;

  //! Number of structures (from config)
  int m_nstructures;
  //! Output freq (from config)
  const int m_outFreq;
  //! acceptance counter
  int m_accepted = 0;
  //! Output model file (from config)
  io::TrajectoryFile m_output;
  //! Output stat file (from config)
  io::StatFile m_outStatsFile;
  //! Switch if move statistics for rigid domains should be logged
  const bool m_logMoveStats;
  //! Inverse temperature (from config)
  const double m_beta;
  //! Restart freq (from config)
  int m_restartFreq;
  //! Restart file - if restart freq >= 0 (from config)
  std::string m_restartFile;

  //! The current sweep loop index
  int m_currentSweepIdx;
  //! session sweep loop counter
  int m_sessionCounter;
  //! current time for dynamic algorithms
  double m_time;

  //! All kernels
  const pairkernels::PairKernelManager &m_kernels;

  //! The stat recorder for the move
  io::MoveStatRecorder m_moveStat;

  ////////////////////////////////////////////////////////////////////////////
  /// Private Membrane Functions
  ////////////////////////////////////////////////////////////////////////////
  virtual int mcSweep() = 0;

  void writeStatHeader() {
    io::writeStatsHeader(m_outStatsFile, "frame", "energy", "accept-ratio",
                         "volume", "energy-connections",
                         m_kernels.getContributionLabels());
  }

  //! Compute the current energy and fill matrix m_energy and set
  //! m_needComputeAll to false.
  void computeAll() {
    m_energy.reset();
    TIMEZONE("computeAll");
    m_interactionComputer->computeAll(m_box, m_forcefield, m_kernels, m_energy);
    m_needComputeAll = false;
  }

  void serializeCore(io::Serializer &serializer) const {
    serializer.append(type(), "type");
    // One by one instead of serializer.append(*m_doms, "m_doms");
    serializer.append(static_cast<size_t>(m_doms->size()), "m_doms->size");
    for (const auto &dom : *m_doms) {
      serializer.append(*dom, "dom");
    }
    serializer.append(m_topologies, "m_topologies");
    serializer.append(m_energy, "m_energy");
    serializer.append(m_box, "m_box");
    // is ref : m_forcefield
    // is ref : m_rng
    serializer.append(*m_interactionComputer, "m_interactionComputer");
    serializer.append(m_needComputeAll, "m_needComputeAll");
    serializer.append(m_nstructures, "m_nstructures");
    serializer.append(m_outFreq, "m_outFreq");
    serializer.append(m_accepted, "m_accepted");
    serializer.append(m_output, "m_output");
    serializer.append(m_outStatsFile, "m_outStatsFile");
    serializer.append(m_logMoveStats, "m_logMoveStats");
    serializer.append(m_beta, "m_beta");
    serializer.append(m_restartFreq, "m_restartFreq");
    serializer.append(m_restartFile, "m_restartFile");
    serializer.append(m_currentSweepIdx, "m_currentSweepIdx");
    serializer.append(m_sessionCounter, "m_sessionCounter");
    serializer.append(m_time, "m_time");
  }

public:
  static std::shared_ptr<domains::Domains>
  RebuildDomains(io::Deserializer &deserializer) {
    const size_t nbDoms = deserializer.restore<size_t>("m_doms->size");
    auto topo = std::make_shared<domains::Domains>(nbDoms);

    for (size_t idx = 0; idx < nbDoms; ++idx) {
      deserializer.access("dom");
      (*topo)[idx] = domains::AbstractDomain::Rebuild(deserializer);
    }

    return topo;
  }

  static std::unique_ptr<AbstractInteractionAlgorithm<double>>
  RebuildInteractionComputer(
      io::Deserializer &deserializer,
      const std::shared_ptr<domains::Domains> &inAllDomains) {
    deserializer.access("m_interactionComputer");
    return AbstractInteractionAlgorithm<double>::Rebuild(deserializer,
                                                         inAllDomains);
  }

  AbstractMcAlgo(io::Deserializer &deserializer,
                 const energy::ForceField &inForcefield,
                 const pairkernels::PairKernelManager &inKernels,
                 util::RNGEngine &inRng)
      : m_doms(RebuildDomains(deserializer)),
        m_topologies(
            deserializer.restore<decltype(m_topologies)>("m_topologies")),
        m_energy(deserializer.restore<decltype(m_energy)>("m_energy")),
        m_box(deserializer.restore<decltype(m_box)>("m_box")),
        m_forcefield(inForcefield), m_rng(inRng),
        m_interactionComputer(RebuildInteractionComputer(deserializer, m_doms)),
        m_needComputeAll(deserializer.restore<decltype(m_needComputeAll)>(
            "m_needComputeAll")),
        m_nstructures(
            deserializer.restore<decltype(m_nstructures)>("m_nstructures")),
        m_outFreq(deserializer.restore<decltype(m_outFreq)>("m_outFreq")),
        m_accepted(deserializer.restore<decltype(m_accepted)>("m_accepted")),
        m_output(deserializer.restore<decltype(m_output)>("m_output")),
        m_outStatsFile(
            deserializer.restore<decltype(m_outStatsFile)>("m_outStatsFile")),
        m_logMoveStats(
            deserializer.restore<decltype(m_logMoveStats)>("m_logMoveStats")),
        m_beta(deserializer.restore<decltype(m_beta)>("m_beta")),
        m_restartFreq(
            deserializer.restore<decltype(m_restartFreq)>("m_restartFreq")),
        m_restartFile(
            deserializer.restore<decltype(m_restartFile)>("m_restartFile")),
        m_currentSweepIdx(deserializer.restore<decltype(m_currentSweepIdx)>(
            "m_currentSweepIdx")),
        m_sessionCounter(deserializer.restore<decltype(m_sessionCounter)>(
            "m_sessionCounter")),
        m_time(deserializer.restore<decltype(m_time)>("m_time")),
        m_kernels(inKernels) {}

  ////////////////////////////////////////////////////////////////////////////
  /// Constructor
  ////////////////////////////////////////////////////////////////////////////

  AbstractMcAlgo(const std::string &inConfigPath, domains::System system,
                 util::RNGEngine &inRng, const util::rvec &inBox,
                 const setup::Config &inConf,
                 const energy::ForceField &inForcefield,
                 const pairkernels::PairKernelManager &inKernels,
                 std::unique_ptr<AbstractInteractionAlgorithm<double>>
                     &&inInteractionComputer)
      : m_doms(system.domains()), m_topologies(system.topologies()),
        m_energy(m_doms->size(), m_doms->size(),
                 inKernels.getNbContributions()),
        m_box(inBox), m_forcefield(inForcefield), m_rng(inRng),
        m_interactionComputer(std::move(inInteractionComputer)),
        m_needComputeAll(true),
        m_nstructures(inConf.value<int>("output.nstructures")),
        m_outFreq(inConf.value<int>("output.freq")),
        m_output(util::appendPathsIfSubRelative(
            inConfigPath, inConf.value<std::string>("output.file"))),
        m_outStatsFile(util::appendPathsIfSubRelative(
            inConfigPath, inConf.value<std::string>("output.stat-file"))),
        // TODO: tremporary HACK to enable log printing of move stats,
        // remove when move stats are in stat file
        m_logMoveStats(
            inConf.experimental_value<bool>("log-move-stats", false)),
        m_beta(constants::natural::refT /
               inConf.value<double>("montecarlo.algorithm-params.temperatur")),
        m_restartFreq(inConf.value<int>("output.restart-freq")),
        m_currentSweepIdx(0), m_sessionCounter(0), m_time(0),
        m_kernels(inKernels) {
    // A restart file must exist it restartFreq >= 0
    if (m_restartFreq >= 0) {
      m_restartFile = util::appendPathsIfSubRelative(
          inConfigPath, inConf.value<std::string>("output.restart-file"));
    } else {
      // No restart
      m_restartFreq = -1;
    }

    if (m_currentSweepIdx == 0) {
      // It is a new simulation, writte headers to files
      computeAll(); // will fill m_energy
      io::writeModel(m_output, (*m_doms), m_box, 0, 0,
                     m_forcefield.beadTypes());
      writeStatHeader();
      io::writeStats(m_outStatsFile, 0, getEnergy(), 0, util::volume(m_box),
                     getEnergyConnections(), getContributions());
      // Real sweep starts at 1
      m_currentSweepIdx = 1;
    }
    m_accepted = 0;
  }

  virtual ~AbstractMcAlgo() {}

  ////////////////////////////////////////////////////////////////////////////
  /// Basic methods
  ////////////////////////////////////////////////////////////////////////////
  int getNbSweeps() const { return m_outFreq * m_nstructures; }

  double getEnergy() const { return m_energy.getTotalEnergy(); }

  double getEnergyConnections() const {
    return m_energy.getTotalEnergyConnections();
  }

  std::vector<double> getContributions() const {
    return m_energy.getTotalContributions();
  }

  const energy::ForceField &getForceField() const { return m_forcefield; }

  double getBeta() const { return m_beta; }

  util::rvec getBox() const { return m_box; }

  virtual double getPressure() const { return 0; }

  int getCurrentSweepIdx() const { return m_currentSweepIdx; }
  int getSessionCounter() const { return m_sessionCounter; }
  void resetSessionCounter() { m_sessionCounter = 0; }

  size_t getNbDomains() const { return m_doms->size(); }

  const io::MoveStatRecorder &getMoveStats() const { return m_moveStat; }

  const pairkernels::PairKernelManager &getKernels() const { return m_kernels; }

  const std::string &getTrajName() const { return m_output.fname(); }
  const std::string &getStatName() const { return m_outStatsFile.fname(); }

  void setTraj(const std::string &fname) {
    m_output = io::TrajectoryFile(fname);
  }
  void setStat(const std::string &fname) {
    m_outStatsFile = io::StatFile(fname);
    writeStatHeader();
  }

  //! Compute the energy for the current topology but the
  //! force field inForceField
  energy::rEnergyMatrix
  computeEnergyForFF(const energy::ForceField &inForceField) const {
    // Trick, we compute if not the same forcefield!
    if (m_forcefield != inForceField) {
      energy::rEnergyMatrix energyForFF(m_doms->size(), m_doms->size(),
                                        m_kernels.getNbContributions());
      m_interactionComputer->computeAll(m_box, inForceField, m_kernels,
                                        energyForFF);
      return energyForFF;
    } else {
      return m_energy;
    }
  }

  //! Exchange coordinate with other.
  //! Here isAccepted will be always true
  //! This will turn m_needComputeAll to true
  //! and cost a computeAll in the call to run method.
  void swapCoordinates(AbstractMcAlgo &other,
                       const std::tuple<bool, double> && /*isAccepted*/) {
    std::swap(m_doms, other.m_doms);
    if (m_energy.sameDimension(other.m_energy)) {
      // Energy matrix follows the domains if possible (to avoid reallocation)
      std::swap(m_energy, other.m_energy);
    } else {
      m_energy = energy::rEnergyMatrix(m_doms->size(), m_doms->size(),
                                       m_kernels.getNbContributions());
      other.m_energy =
          energy::rEnergyMatrix(other.m_doms->size(), other.m_doms->size(),
                                other.m_kernels.getNbContributions());
    }
    // The interaction computer too
    std::swap(m_interactionComputer, other.m_interactionComputer);
    std::swap(m_box, other.m_box);

    m_needComputeAll = true;
    other.m_needComputeAll = true;
  }

  //! Exchange coordinate with other.
  //! Use the energy matrices from isAcceptedWithEnergy
  //! as current state.
  //! isAcceptedWithEnergy contains
  //! - true
  //! - other.forcefield(this)
  //! - this.forcefield(other)
  void swapCoordinates(
      AbstractMcAlgo &other,
      std::tuple<bool, double, energy::rEnergyMatrix, energy::rEnergyMatrix>
          &&isAcceptedWithEnergy) {
    std::swap(m_doms, other.m_doms);
    // Use the matrices that has been computed in the accept function
    m_energy = std::move(std::get<3>(isAcceptedWithEnergy));
    other.m_energy = std::move(std::get<2>(isAcceptedWithEnergy));
    // The interaction computer too
    std::swap(m_interactionComputer, other.m_interactionComputer);
    std::swap(m_box, other.m_box);
  }

  void setNStructures(int nStruc) { m_nstructures = nStruc; }

  //! Must be provided by sub class
  virtual double accept(const double deltaE) = 0;

  ////////////////////////////////////////////////////////////////////////////
  /// Run
  ////////////////////////////////////////////////////////////////////////////

  //! Note that the sweep 0 is skipped.
  //! The first sweep loop index is 1
  void run(const int nbSweepsToProceed) {
    TIMEZONE("mcAlgo");

    // Recompute the energy matrix if needed
    if (m_needComputeAll) {
      computeAll();
    }
    const int lastSweepIdx = m_currentSweepIdx + nbSweepsToProceed;
    for (; m_currentSweepIdx < lastSweepIdx;
         ++m_currentSweepIdx, ++m_sessionCounter) {
      // first write restart file if required. This ensures that there is at
      // least one overlap frame when we restart a trajectory to check if it
      // indeed continuous the old one.
      if (m_restartFreq != -1 &&
          ((m_currentSweepIdx + 1) % m_restartFreq == 0)) {
        util::Log("writing restart file at sweep {}\n", m_currentSweepIdx);
        io::writeRestart(m_restartFile, *this, m_rng);
        m_output.flush();
        m_outStatsFile.flush();
      }

      TIMEZONE("SweepLoop");
      if (m_logMoveStats) {
        util::Log("MOVE-TYPE-STATS: sweep: {}\n", m_currentSweepIdx);
      }
      // accepted volumeMoves will also be returned as 1 here.
      m_accepted += mcSweep();
      // Write output if needed
      if ((m_currentSweepIdx % m_outFreq == 0)) {
        io::writeModel(m_output, (*m_doms), m_box,
                       m_currentSweepIdx / m_outFreq, m_time,
                       m_forcefield.beadTypes());
        io::writeStats(
            m_outStatsFile, m_currentSweepIdx / m_outFreq,
            m_energy.getTotalEnergy(),
            m_accepted / static_cast<double>((m_outFreq)*m_doms->size()),
            util::volume(m_box), m_energy.getTotalEnergyConnections(),
            m_energy.getTotalContributions());
        m_accepted = 0;
      }
    }
  }

  virtual std::string type() const = 0;
};

//! This class provide a method to the
//! template method pattern using a simple function.
//! This class allows to use already define simple funtions
//! under the AbstractMcAlgo class.
template <class ParentMcAlgo, double Acceptor(double, double)>
class McAlgoWithAcceptFunc : public ParentMcAlgo {
public:
  using ParentMcAlgo::ParentMcAlgo;

  double accept(const double deltaE) override {
    return Acceptor(deltaE, ParentMcAlgo::m_beta);
  }

  std::string type() const final { return Type(); }

  static std::string Type() { return __PRETTY_FUNCTION__; }
};

} // namespace mc

#endif
