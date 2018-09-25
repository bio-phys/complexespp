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
#include "domains/rigiddomain.h"
#include "energy/forcefield.h"
#include "mc/accept_func.h"
#include "setup/config.h"
#include "util/array.h"
#include "util/random.h"

#include "mc/nvt.h"

using util::RNGEngine;

namespace mc {

int NVTMC::mcSweep() {
  TIMEZONE("mcSweep");
  auto movement_counter = 0;
  auto trans_moves = 0;
  auto trans_accepted = 0;
  auto rot_moves = 0;
  auto rot_accepted = 0;
  auto accepted = 0;
  auto uniform =
      std::uniform_int_distribution<>{0, static_cast<int>(m_doms->size()) - 1};

  for (auto i = 0u; i < m_doms->size(); ++i) {
    movement_counter += 1;
    auto index = uniform(m_rng);
    domains::MovedDomain newState;
    // Move it
    {
      TIMEZONE("move");
      newState = m_doms->at(index)->move(m_box, m_rng);
    }
    if (newState.moveSucceed()) {
      //      util::Log("MOVE-TYPE-STATS {}\n", newState.type());
      if (newState.type() == "trans") {
        trans_moves += 1;
      } else if (newState.type() == "rot") {
        rot_moves += 1;
      }
      // Keep ownership of the unmoved state
      std::unique_ptr<domains::AbstractDomain> previousState =
          std::move(m_doms->at(index));
      m_doms->at(index) = newState.releaseDomain();

      // Inform the interaction-algorithm about the uptdate to make
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
        // Use access over the rows (since m_energy is likely to be RowMajor)
        enDelta += enNew.getEnergy(0, j) - m_energy.getEnergy(j, index);
      }

      if (std::isnan(enDelta) || std::isinf(enDelta)) {
        util::Log("Invalid energy ({}), this might come from a collision.\n",
                  enDelta);
        util::Log("The programm will continue but is no longer valid.\n");
      }

      // Check if the move is accepted
      if (std::uniform_real_distribution<>{0, 1}(m_rng) < accept(enDelta)) {
        // Count accepted times
        if (newState.type() == "trans") {
          trans_accepted += 1;
        } else if (newState.type() == "rot") {
          rot_accepted += 1;
        }
        ++accepted;
        // Update the energies matrix with the newest values
        m_energy.replaceRowAndCol(index, enNew);
        m_moveStat.addEvent(index, m_doms->at(index)->name(), true, true);
      } else {
        // Not accepted, rollback to the unmoved state
        m_doms->at(index) = std::move(previousState);
        // Inform the inter-algo about the update to make
        TIMEZONE("updateDomain");
        m_interactionComputer->updateDomain(index);
        m_moveStat.addEvent(index, m_doms->at(index)->name(), true, false);
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
  if (m_logMoveStats) {
    util::Log("MOVE-TYPE-STATS: trans_moves: {}; rot_moves {}; total: {}\n",
              trans_moves, rot_moves, movement_counter);
    util::Log("MOVE-TYPE-STATS: trans_accepted: {}; rot_accepted {}\n",
              trans_accepted, rot_accepted);
  }

  for (auto &topology : m_topologies) {
    if (topology.move(m_box, m_doms.get(), &m_rng)) {
      // We inform the container that these domains have moved
      for (auto domainId : topology.domainIds()) {
        m_interactionComputer->updateDomain(domainId);
      }
    }
  }

  return accepted;
}

//! This function create a McAlgo (nvt) based on the
//! accept function specified in the parameter/config.
std::unique_ptr<NVTMC>
McNVTBuild(const std::string &inConfigPath, domains::System system,
           RNGEngine &rng, const util::rvec &box, const setup::Config &conf,
           const energy::ForceField &forcefield,
           const pairkernels::PairKernelManager &inKernels,
           std::unique_ptr<AbstractInteractionAlgorithm<double>>
               &&interactionComputer) {
  if (conf.value<std::string>("montecarlo.algorithm-params.accept-func") ==
      "metropolis") {
    return std::make_unique<McAlgoWithAcceptFunc<NVTMC, metropolisAccept>>(
        inConfigPath, system, rng, box, conf, forcefield, inKernels,
        std::move(interactionComputer));
  } else if (conf.value<std::string>(
                 "montecarlo.algorithm-params.accept-func") == "glauber") {
    return std::make_unique<McAlgoWithAcceptFunc<NVTMC, glauberAccept>>(
        inConfigPath, system, rng, box, conf, forcefield, inKernels,
        std::move(interactionComputer));
  } else if (conf.value<std::string>(
                 "montecarlo.algorithm-params.accept-func") == "dynamic") {
    return std::make_unique<McAlgoWithAcceptFunc<NVTMC, dynamicAccept>>(
        inConfigPath, system, rng, box, conf, forcefield, inKernels,
        std::move(interactionComputer));
  } else if (conf.value<std::string>(
                 "montecarlo.algorithm-params.accept-func") == "always") {
    return std::make_unique<McAlgoWithAcceptFunc<NVTMC, alwaysAccept>>(
        inConfigPath, system, rng, box, conf, forcefield, inKernels,
        std::move(interactionComputer));
  } else {
    throw std::invalid_argument(fmt::format(
        "{} -- Unknown acceptance function",
        conf.value<std::string>("montecarlo.algorithm.accept-func")));
  }
  // return std::unique_ptr<McAlgo>();
}
} // namespace mc
