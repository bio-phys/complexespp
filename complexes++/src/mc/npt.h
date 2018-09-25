// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#ifndef MC_NPT_H
#define MC_NPT_H

#include "abstractmcalgo.h"
#include "domains/abstractdomain.h"
#include "mc/accept_func.h"

class ForceField;
class Config;

namespace mc {

class NPTMC : public AbstractMcAlgo {
 protected:
  double getPressure() const final { return m_pressure; }
  int mcSweep() final;
  int volumeMove();

 public:
  NPTMC(const std::string& inConfigPath, domains::System system,
        util::RNGEngine& inRng, const util::rvec& inBox,
        const setup::Config& inConf, const energy::ForceField& inForcefield,
        const pairkernels::PairKernelManager& inKernels,
        std::unique_ptr<AbstractInteractionAlgorithm<double>>&&
            inInteractionComputer);

  virtual ~NPTMC() {}

  void serialize(io::Serializer& serializer) const final {
    AbstractMcAlgo::serializeCore(serializer);
    serializer.append(m_pressure, "m_pressure");
    serializer.append(m_dV, "m_dV");
    serializer.append(m_V, "m_V");
  }

  NPTMC(io::Deserializer& deserializer, const energy::ForceField& inForcefield,
        const pairkernels::PairKernelManager& inKernels, util::RNGEngine& inRng)
      : AbstractMcAlgo(deserializer, inForcefield, inKernels, inRng),
        m_pressure(deserializer.restore<decltype(m_pressure)>("m_pressure")),
        m_dV(deserializer.restore<decltype(m_dV)>("m_dV")),
        m_V(deserializer.restore<decltype(m_V)>("m_V")) {}

 private:
  double m_pressure;  // should be in kt / \AA^3
  double m_dV;
  double m_V;
  bool m_verbose;
};
using NPTMC_metropolisAccept =
    mc::McAlgoWithAcceptFunc<NPTMC, mc::metropolisAccept>;
REBUILDER_REGISTER(NPTMC_metropolisAccept);
using NPTMC_glauberAccept = mc::McAlgoWithAcceptFunc<NPTMC, mc::glauberAccept>;
REBUILDER_REGISTER(NPTMC_glauberAccept);
using NPTMC_dynamicAccept = mc::McAlgoWithAcceptFunc<NPTMC, mc::dynamicAccept>;
REBUILDER_REGISTER(NPTMC_dynamicAccept);
using NPTMC_alwaysAccept = mc::McAlgoWithAcceptFunc<NPTMC, mc::alwaysAccept>;
REBUILDER_REGISTER(NPTMC_alwaysAccept);

std::unique_ptr<NPTMC> McNPTBuild(
    const std::string& inConfigPath, domains::System system,
    util::RNGEngine& rng, const util::rvec& box, const setup::Config& conf,
    const energy::ForceField& forcefield,
    const pairkernels::PairKernelManager& inKernels,
    std::unique_ptr<AbstractInteractionAlgorithm<double>>&&
        interactionComputer);
}  // namespace mc

#endif  // MC_NPT_H
