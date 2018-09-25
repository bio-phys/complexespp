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
#include "domains/rigiddomain.h"
#include "io/ffparameters.h"
#include "pairkernels/ljhpairkernel.h"
#include "gtest/gtest.h"

domains::Rigid singleAtomRigid(const int id, const domains::Bead bead,
                               const double charge, const util::rvec &pos) {
  const auto Dtrans = util::rvec(.2, .2, .2);
  const auto phi = 2.0;
  auto dom = domains::Rigid(
      "none", 0, id, std::vector<domains::Bead>(1, domains::Bead(bead)),
      std::vector<double>{charge},
      std::vector<domains::BeadChainID>{domains::BeadChainID("A", 1)},
      domains::Connections(0), Dtrans, phi);

  auto xyz = util::rArray(1, 3);
  for (auto i = 0; i < xyz.cols(); ++i) {
    xyz(0, i) = pos[i];
  }
  dom.setXyz(std::move(xyz));
  return dom;
}

class RigidTest : public testing::Test {
public:
  RigidTest() {}
  virtual ~RigidTest() {}

protected:
  virtual void SetUp() {
    const auto beadTypes = io::readBeadTypes(dataDir + "bead-types");
    // r0 is 2^(1/6) minimum of energy function
    auto third_of_r0 = std::sqrt(std::pow(2, 1. / 3.) / 3);
    dom1 = singleAtomRigid(1, domains::findBeadID("ALA", beadTypes), 0,
                           util::rvec(0, 0, 0));
    dom2 = singleAtomRigid(2, domains::findBeadID("ALA", beadTypes), 0,
                           util::rvec(third_of_r0, third_of_r0, third_of_r0));
  }
  virtual void TearDown() {}

  double calc_energy(const domains::Rigid &a, const domains::Rigid &b) {
    if (a.id() != b.id()) {
      energy::rEnergyMatrix array(
          1, 1, pairkernels::LJHPairKernel::LJHNbContributions);
      omp_lock_t energyMatrixMutex;
      omp_init_lock(&energyMatrixMutex);
      {
        energy::EnergyMatrixBuffer<> bufferEnergyArray(array,
                                                       energyMatrixMutex);
        pairkernels::LJHPairKernel().AbstractPairKernel::compute(
            bufferEnergyArray.getAccesser(0, 0), a, b, box, forcefield);
      }
      omp_destroy_lock(&energyMatrixMutex);
      return array.getEnergy(0, 0);
    } else {
      return 0;
    }
  }

  domains::Rigid dom1 =
      singleAtomRigid(1, domains::Bead(0), 0, util::rvec(0, 0, 0));
  domains::Rigid dom2 =
      singleAtomRigid(2, domains::Bead(0), 0, util::rvec(0, 0, 0));

private:
  const energy::ForceField forcefield = dummy_forcefield(1, -1, 1, 1, 1);
  util::rvec box = util::rvec(10, 10, 10);
};

TEST_F(RigidTest, energy_symmetric) {
  auto en = calc_energy(dom1, dom2);
  auto en2 = calc_energy(dom2, dom1);
  EXPECT_DOUBLE_EQ(en, en2);
  EXPECT_DOUBLE_EQ(-1, en);
}

TEST_F(RigidTest, energySelf) { EXPECT_DOUBLE_EQ(0, calc_energy(dom1, dom1)); }
