// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#include "mc/mc.h"
#include "mc/accept_func.h"
#include "util/random.h"

#include "mc/abstractmcalgo.h"
#include "mc/npt.h"
#include "mc/nvt.h"

namespace mc {

//! This function create a McAlgo based on the given config (conf)
//! It could return a Brownian or normal MC.
//! It also check for collision using parameters to know if
//! advanced method should be used.
std::unique_ptr<AbstractMcAlgo> McBuild(
    const std::string& inConfigPath, domains::System system,
    util::RNGEngine& rng, const util::rvec& box, const setup::Config& conf,
    const energy::ForceField& forcefield,
    const pairkernels::PairKernelManager& inKernels) {
  // to ensure that no fatal overlaps are present when we calculate the energy.
  // This can happen with gaussian domains for which it is acceptable to all
  // overlap in the initial state.
  std::unique_ptr<AbstractInteractionAlgorithm<double>> interactionComputer =
      buildInteractionComputer<double>(system.domains(), box, conf);

  auto const& algorithm = conf.value<std::string>("montecarlo.algorithm");
  // 'normal' is still allowed for backwards compatibility
  if (algorithm == "nvt" || algorithm == "normal") {
    return McNVTBuild(inConfigPath, system, rng, box, conf, forcefield,
                      inKernels, std::move(interactionComputer));
  } else if (algorithm == "npt") {
    return McNPTBuild(inConfigPath, system, rng, box, conf, forcefield,
                      inKernels, std::move(interactionComputer));
  } else {
    throw std::invalid_argument(
        fmt::format("{} -- unkown montecarlo algorithm", algorithm));
  }
}

void McRerun(const std::string& baseDir, std::shared_ptr<domains::Domains> doms,
             const util::rvec& box, const setup::Config& conf,
             const energy::ForceField& forcefield,
             const pairkernels::PairKernelManager& inKernels) {
  auto reader = io::readTrajectory(
      doms, box,
      util::appendPaths(baseDir, conf.value<std::string>("output.file")));
  io::StatFile stat(
      util::appendPaths(baseDir, conf.value<std::string>("output.stat-file")));
  io::writeStatsHeader(stat, "frame", "energy", "energy-connections",
                       inKernels.getContributionLabels());

  std::unique_ptr<AbstractInteractionAlgorithm<double>> interactionComputer =
      buildInteractionComputer<double>(doms, box, conf);

  auto i = 0;
  while (reader->hasNextFrame()) {
    // Go to next frame does not change domains shared ptr
    // 'doms' is still valid for all loops
    reader->nextFrame();
    interactionComputer->resetAllDomains(box);
    auto energy = energy::rEnergyMatrix(doms->size(), doms->size(),
                                        inKernels.getNbContributions());
    interactionComputer->computeAll(box, forcefield, inKernels, energy);
    io::writeStats(stat, i++, energy.getTotalEnergy(),
                   energy.getTotalEnergyConnections(),
                   energy.getTotalContributions());
  }
}
}  // namespace mc
