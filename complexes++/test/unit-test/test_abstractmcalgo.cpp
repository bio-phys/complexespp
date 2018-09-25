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
#include "complexes_test.h"
#include "gtest/gtest.h"

#include "domains/abstractdomain.h"
#include "io/cplx.h"
#include "mc/abstractmcalgo.h"
#include "mc/accept_func.h"
#include "mc/nvt.h"
#include "setup/config.h"

class TestMC : public mc::AbstractMcAlgo {
protected:
  int mcSweep() final {
    // ensure getEnergy returns the same energy as stored in the energy matrix
    EXPECT_DOUBLE_EQ(getEnergy(), m_energy.getTotalEnergy());

    for (auto i = 0u; i < m_doms->size(); ++i) {
      auto movedDomain = m_doms->at(i)->move(m_box, m_rng);
      if (movedDomain.moveSucceed()) {
        m_doms->at(i) = movedDomain.releaseDomain();
        // Inform the interaction-algorithm about the uptdate to make
        m_interactionComputer->updateDomain(i);

        // Compute the new energy with this move
        energy::rEnergyMatrix enNew(1UL, m_doms->size(),
                                    getKernels().getNbContributions());
        {
          TIMEZONE("computeForOneDomain");
          m_interactionComputer->computeForOneDomain(i, m_box, m_forcefield,
                                                     getKernels(), enNew);
        }

        // update energy matrix
        m_energy.replaceRowAndCol(i, enNew);
      }
    }

    return static_cast<int>(m_doms->size());
  }

public:
  using AbstractMcAlgo::AbstractMcAlgo;

  void serialize(io::Serializer &serializer) const final {
    AbstractMcAlgo::serializeCore(serializer);
  }

  virtual ~TestMC() {}
};

std::unique_ptr<TestMC>
McTestBuild(std::shared_ptr<domains::Domains> doms, util::RNGEngine &rng,
            const util::rvec &box, const setup::Config &conf,
            const energy::ForceField &forcefield,
            pairkernels::PairKernelManager &inKernels,
            std::unique_ptr<AbstractInteractionAlgorithm<double>>
                &&interactionComputer) {
  return std::make_unique<
      mc::McAlgoWithAcceptFunc<TestMC, mc::metropolisAccept>>(
      "", domains::System(doms, std::vector<domains::Topology>()), rng, box,
      conf, forcefield, inKernels, std::move(interactionComputer));
}

class AbstractMcAlgoTest : public testing::Test {
public:
  AbstractMcAlgoTest() {}
  virtual ~AbstractMcAlgoTest() {}

  std::shared_ptr<domains::System> m_doms;
  std::unique_ptr<AbstractInteractionAlgorithm<double>> m_interactionComputer;
  std::unique_ptr<TestMC> m_algo;

  pairkernels::PairKernelManager m_kernels;

protected:
  virtual void SetUp() {
    // Carefull in future one might want to have reproducible random numbers
    std::random_device rd;
    m_rng = util::RNGEngine(rd());

    m_doms = std::make_shared<domains::System>(
        io::readCPLX(dataDir + "threeDomain.cplx", m_forcefield.beadTypes()));

    m_kernels =
        io::readKernelMapping(dataDir + "threeDomain.cplx", *m_doms->domains());

    m_algo = McTestBuild(
        m_doms->domains(), m_rng, m_box, m_config, m_forcefield, m_kernels,
        buildInteractionComputer<double>(m_doms->domains(), m_box, m_config));
  }
  virtual void TearDown() {}

private:
  const energy::ForceField m_forcefield = dummy_forcefield(1, -1, 1, 1, 1);
  const util::rvec m_box = util::rvec(10, 10, 10);
  util::RNGEngine m_rng;
  setup::Config m_config = setup::Config(dataDir + "abstractalgo.conf");
};

TEST_F(AbstractMcAlgoTest, run_10_sweeps) {
  // 0 is reserved for the initial state
  EXPECT_EQ(1, m_algo->getCurrentSweepIdx());
  m_algo->run(10);
  EXPECT_EQ(11, m_algo->getCurrentSweepIdx());
  EXPECT_DOUBLE_EQ(1., m_algo->getBeta());
  EXPECT_EQ(10, m_algo->getSessionCounter());
  m_algo->resetSessionCounter();
  EXPECT_EQ(0, m_algo->getSessionCounter());
}
