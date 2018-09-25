// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#include "mc/npt.h"
#include "domains/rigiddomain.h"
#include "energy/forcefield.h"
#include "mc/accept_func.h"
#include "mc/nvt.h"
#include "setup/config.h"
#include "util/array.h"
#include "util/random.h"

using util::RNGEngine;

namespace mc {

NPTMC::NPTMC(const std::string& inConfigPath, domains::System system,
             util::RNGEngine& inRng, const util::rvec& inBox,
             const setup::Config& inConf,
             const energy::ForceField& inForcefield,
             const pairkernels::PairKernelManager& inKernels,
             std::unique_ptr<AbstractInteractionAlgorithm<double>>&&
                 inInteractionComputer)
    : AbstractMcAlgo(inConfigPath, system, inRng, inBox, inConf, inForcefield,
                     inKernels, std::move(inInteractionComputer)),
      m_pressure(constants::conversions::barTOktperCubAA *
                 inConf.value<double>("montecarlo.algorithm-params.pressure")),
      m_dV(inConf.value<double>("montecarlo.algorithm-params.dV")),
      m_V(util::volume(inBox)),
      m_verbose(inConf.value<bool>("montecarlo.algorithm-params.verbose")) {}

int NPTMC::volumeMove() {
  const auto oldEN = getEnergy();
  const auto oldBox = m_box;
  // 1. calculate new volume
  const auto dV = std::uniform_real_distribution<>{-m_dV, m_dV}(m_rng);
  m_V = util::volume(m_box);
  const auto newV = m_V + dV;

  if (newV < 0) {  // discard move without further evaluation
    util::Log(
        "New Volume would be negative. dV is in the same magnitude as V "
        "(or bigger)\n"
        "Volume-move will be discarded.\n");
    const auto notAccepted = false;
    return notAccepted;
  }
  // 2. determine scaling factor
  const auto scale = cbrt(newV / m_V);
  // 3. scale box and coordinates
  const auto new_box = util::scale(m_box, scale);
  auto oldXyz = std::vector<util::rArray>();
  oldXyz.reserve(m_doms->size());
  for (auto& dom : *m_doms) {
    oldXyz.push_back(dom->xyz());
    auto _xyz = dom->xyz();
    dom->setXyz(util::scaleCentroid(_xyz, scale));
  }
  m_interactionComputer->resetAllDomains(new_box);

  // 6. evaluate total energy
  // use new box for energy calculations!
  m_box = new_box;
  computeAll();

  // 6. accept reject, From Daan Frenkel (2002) eq 5.4.11, Understanding
  // Molecular Simulations
  const auto N = static_cast<int>(m_doms->size());
  const auto newEN = getEnergy();
  const auto deltaE = newEN - oldEN;
  const auto acceptProb =
      accept(deltaE + m_pressure * dV - N / m_beta * std::log(newV / m_V));

  const auto isAccepted =
      std::uniform_real_distribution<>{0, 1}(m_rng) < acceptProb;
  if (isAccepted) {
    m_V += dV;
    // the rest was already updated to calculate the energies
  } else {
    for (auto& dom : *m_doms) {
      dom->setXyz(oldXyz[dom->id()]);
    }
    m_box = oldBox;
    m_interactionComputer->resetAllDomains(m_box);
    computeAll();
  }

  return isAccepted;
}

int NPTMC::mcSweep() {
  TIMEZONE("mcSweep");
  auto movement_counter = 0;
  auto trans_moves = 0;
  auto trans_accepted = 0;
  auto rot_moves = 0;
  auto rot_accepted = 0;
  auto volume_move_counter = 0;
  auto volume_move_acceptance_counter = 0;
  auto accepted = 0;
  // TODO make input somewhere
  auto volume_move_frq = 1;
  auto uniform =
      std::uniform_int_distribution<>{0, static_cast<int>(m_doms->size()) - 1};

  for (auto i = 0u; i < (m_doms->size() + volume_move_frq); i++) {
    if (std::uniform_real_distribution<>{0, 1}(m_rng) <
        1.0 / static_cast<double>(m_doms->size() + volume_move_frq)) {
      volume_move_counter += 1;
      auto volume_move_accepted = volumeMove();
      if (volume_move_accepted) {
        volume_move_acceptance_counter += 1;
      }
    } else {
      movement_counter += 1;
      auto index = uniform(m_rng);
      domains::MovedDomain newState;
      // Move it
      {
        TIMEZONE("move");
        newState = m_doms->at(index)->move(m_box, m_rng);
      }
      // const auto& type = newState.type();
      if (newState.moveSucceed()) {
        if (newState.type() == "trans") {
          trans_moves += 1;
        } else if (newState.type() == "rot") {
          rot_moves += 1;
        }
        std::unique_ptr<domains::AbstractDomain> previousState =
            std::move(m_doms->at(index));
        m_doms->at(index) = newState.releaseDomain();
        // Inform the interaction-algorithm about the update to make
        m_interactionComputer->updateDomain(index);
        // Compute the new energy with this move
        energy::rEnergyMatrix enNew(1UL, m_doms->size(),
                                    m_kernels.getNbContributions());
        {
          TIMEZONE("computeForOneDomain");
          m_interactionComputer->computeForOneDomain(index, m_box, m_forcefield,
                                                     m_kernels, enNew);
        }
        // Compute the delta from previous energies and new ones
        auto enDelta = 0.0;
        for (auto j = 0u; j < m_doms->size(); ++j) {
          // Use access over the rows (since m_energy is likely to be
          // RowMajor)
          enDelta += enNew.getEnergy(0, j) - m_energy.getEnergy(j, index);
        }

        if (std::isnan(enDelta) || std::isinf(enDelta)) {
          util::Log("Invalid energy ({}), this might come from a collision.\n",
                    enDelta);
          util::Log("The programm will continue but is no longer valid.\n");
        }
        // Check if the move is accepted
        if (std::uniform_real_distribution<>{0, 1}(m_rng) < accept(enDelta)) {
          if (newState.type() == "trans") {
            trans_accepted += 1;
          } else if (newState.type() == "rot") {
            rot_accepted += 1;
          }
          // Count accepted times
          ++accepted;
          // Update the energies matrix with the newest values
          m_energy.replaceRowAndCol(index, enNew);
          m_moveStat.addEvent(index, m_doms->at(index)->name(), false, false);
        } else {
          // Not accepted, rollback to the unmoved state
          m_doms->at(index) = std::move(previousState);
          // Inform the inter-algo about the update to make
          TIMEZONE("updateDomain");
          m_interactionComputer->updateDomain(index);
        }
      } else {
        // Check if the move is accepted for a delta energy = 0
        if (std::uniform_real_distribution<>{0, 1}(m_rng) < accept(0.)) {
          // Count accepted times, but nothing to do
          ++accepted;
        }
        m_moveStat.addEvent(index, m_doms->at(index)->name(), false, false);
      }
    }
    // update acceptance counters
  }
  // Log acceptance counter
  if (m_verbose) {
    util::Log("MOVE-TYPE-STATS: trans_moves: {}; rot_moves {}; total: {}\n",
              trans_moves, rot_moves, movement_counter);
    util::Log("MOVE-TYPE-STATS: trans_accepted: {}; rot_accepted {}\n",
              trans_accepted, rot_accepted);
    util::Log("MOVE-TYPE-STATS: volume_moves: {}; volumes_moves_accepted {}\n",
              volume_move_counter, volume_move_acceptance_counter);
  }
  return accepted;
}
//! This function create a McAlgo (normal) based on the
//! accept function specified in the parameter/config.
std::unique_ptr<NPTMC> McNPTBuild(
    const std::string& inConfigPath, domains::System system,
    util::RNGEngine& rng, const util::rvec& box, const setup::Config& conf,
    const energy::ForceField& forcefield,
    const pairkernels::PairKernelManager& inKernels,
    std::unique_ptr<AbstractInteractionAlgorithm<double>>&&
        interactionComputer) {
  if (conf.value<std::string>("montecarlo.algorithm-params.accept-func") ==
      "metropolis") {
    return std::make_unique<McAlgoWithAcceptFunc<NPTMC, metropolisAccept>>(
        inConfigPath, system, rng, box, conf, forcefield, inKernels,
        std::move(interactionComputer));
  } else if (conf.value<std::string>(
                 "montecarlo.algorithm-params.accept-func") == "glauber") {
    return std::make_unique<McAlgoWithAcceptFunc<NPTMC, glauberAccept>>(
        inConfigPath, system, rng, box, conf, forcefield, inKernels,
        std::move(interactionComputer));
  } else if (conf.value<std::string>(
                 "montecarlo.algorithm-params.accept-func") == "dynamic") {
    return std::make_unique<McAlgoWithAcceptFunc<NPTMC, dynamicAccept>>(
        inConfigPath, system, rng, box, conf, forcefield, inKernels,
        std::move(interactionComputer));
  } else if (conf.value<std::string>(
                 "montecarlo.algorithm-params.accept-func") == "always") {
    return std::make_unique<McAlgoWithAcceptFunc<NPTMC, alwaysAccept>>(
        inConfigPath, system, rng, box, conf, forcefield, inKernels,
        std::move(interactionComputer));
  } else {
    throw std::invalid_argument(fmt::format(
        "{} -- Unknown acceptance function",
        conf.value<std::string>("montecarlo.algorithm.accept-func")));
  }
}
}  // namespace mc
