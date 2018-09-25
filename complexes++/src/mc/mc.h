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
std::unique_ptr<AbstractMcAlgo>
McBuild(const std::string &inConfigPath, domains::System system,
        util::RNGEngine &rng, const util::rvec &box, const setup::Config &conf,
        const energy::ForceField &forcefield,
        const pairkernels::PairKernelManager &inKernels);

void McRerun(const std::string &baseDir, std::shared_ptr<domains::Domains> doms,
             const util::rvec &box, const setup::Config &conf,
             const energy::ForceField &forcefield,
             const pairkernels::PairKernelManager &inKernels);

} // namespace mc

#endif // MC_H
