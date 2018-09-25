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
  NPTMC(const std::string &inConfigPath, domains::System system,
        util::RNGEngine &inRng, const util::rvec &inBox,
        const setup::Config &inConf, const energy::ForceField &inForcefield,
        const pairkernels::PairKernelManager &inKernels,
        std::unique_ptr<AbstractInteractionAlgorithm<double>>
            &&inInteractionComputer);

  virtual ~NPTMC() {}

  void serialize(io::Serializer &serializer) const final {
    AbstractMcAlgo::serializeCore(serializer);
    serializer.append(m_pressure, "m_pressure");
    serializer.append(m_dV, "m_dV");
    serializer.append(m_V, "m_V");
  }

  NPTMC(io::Deserializer &deserializer, const energy::ForceField &inForcefield,
        const pairkernels::PairKernelManager &inKernels, util::RNGEngine &inRng)
      : AbstractMcAlgo(deserializer, inForcefield, inKernels, inRng),
        m_pressure(deserializer.restore<decltype(m_pressure)>("m_pressure")),
        m_dV(deserializer.restore<decltype(m_dV)>("m_dV")),
        m_V(deserializer.restore<decltype(m_V)>("m_V")) {}

private:
  double m_pressure; // should be in kt / \AA^3
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

std::unique_ptr<NPTMC>
McNPTBuild(const std::string &inConfigPath, domains::System system,
           util::RNGEngine &rng, const util::rvec &box,
           const setup::Config &conf, const energy::ForceField &forcefield,
           const pairkernels::PairKernelManager &inKernels,
           std::unique_ptr<AbstractInteractionAlgorithm<double>>
               &&interactionComputer);
} // namespace mc

#endif // MC_NPT_H
