// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
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

  void serialize(io::Serializer& serializer) const final {
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

std::unique_ptr<NVTMC> McNVTBuild(
    const std::string& inConfigPath, domains::System system,
    util::RNGEngine& rng, const util::rvec& box, const setup::Config& conf,
    const energy::ForceField& forcefield,
    const pairkernels::PairKernelManager& inKernels,
    std::unique_ptr<AbstractInteractionAlgorithm<double>>&&
        interactionComputer);
}  // namespace mc

#endif  // MC_NVT_H
