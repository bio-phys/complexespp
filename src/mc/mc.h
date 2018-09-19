// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#ifndef MC_H
#define MC_H

#include "domains/abstractdomain.h"
#include "energy/forcefield.h"
#include "interactions/interactionbuilder.h"
#include "mc/abstractmcalgo.h"
#include "pairkernels/pairkernelmanager.h"
#include "setup/config.h"
#include "util/random.h"

namespace mc {
std::unique_ptr<AbstractMcAlgo> McBuild(
    const std::string& inConfigPath, domains::System system,
    util::RNGEngine& rng, const util::rvec& box, const setup::Config& conf,
    const energy::ForceField& forcefield,
    const pairkernels::PairKernelManager& inKernels);

void McRerun(const std::string& baseDir, std::shared_ptr<domains::Domains> doms,
             const util::rvec& box, const setup::Config& conf,
             const energy::ForceField& forcefield,
             const pairkernels::PairKernelManager& inKernels);

}  // namespace mc

#endif  // MC_H
