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
#ifndef MC_NVT_H
#define MC_NVT_H

#include "domains/abstractdomain.h"
#include "domains/moveddomain.h"
#include "mc/abstractmcalgo.h"

class RNGEnging;
class ForceField;
class Config;

namespace mc {

class NVTMC : public AbstractMcAlgo {
protected:
  ////////////////////////////////////////////////////////////////////////////
  /// Sweep
  ////////////////////////////////////////////////////////////////////////////

  int mcSweep() final;

public:
  using AbstractMcAlgo::AbstractMcAlgo;

  void serialize(io::Serializer &serializer) const final {
    AbstractMcAlgo::serializeCore(serializer);
  }

  virtual ~NVTMC() {}
};
using NVTMC_metropolisAccept =
    mc::McAlgoWithAcceptFunc<NVTMC, mc::metropolisAccept>;
REBUILDER_REGISTER(NVTMC_metropolisAccept);
using NVTMC_glauberAccept = mc::McAlgoWithAcceptFunc<NVTMC, mc::glauberAccept>;
REBUILDER_REGISTER(NVTMC_glauberAccept);
using NVTMC_dynamicAccept = mc::McAlgoWithAcceptFunc<NVTMC, mc::dynamicAccept>;
REBUILDER_REGISTER(NVTMC_dynamicAccept);
using NVTMC_alwaysAccept = mc::McAlgoWithAcceptFunc<NVTMC, mc::alwaysAccept>;
REBUILDER_REGISTER(NVTMC_alwaysAccept);

std::unique_ptr<NVTMC>
McNVTBuild(const std::string &inConfigPath, domains::System system,
           util::RNGEngine &rng, const util::rvec &box,
           const setup::Config &conf, const energy::ForceField &forcefield,
           const pairkernels::PairKernelManager &inKernels,
           std::unique_ptr<AbstractInteractionAlgorithm<double>>
               &&interactionComputer);
} // namespace mc

#endif // MC_NVT_H
