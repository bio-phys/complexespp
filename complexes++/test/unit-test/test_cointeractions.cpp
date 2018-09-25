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
#include "gtest/gtest.h"

#include "cutoffgrid/cosparsegridcontainer.h"
#include "cutoffgrid/cutoffinteractions.h"
#include "domains/abstractdomain.h"
#include "domains/rigiddomain.h"
#include "energy/energy.h"
#include "interactions/fullinteractions.h"
#include "util/moves.h"
#include "util/util.h"

#include "complexes_test.h"

#include "cotestingdomain.h"

TEST(CUTOFF_TEST, testInteractionAlgorithms) {
  const double coRadius = 10;
  const auto boxSize = util::rvec(100, 200, 300);
  const int nbDomains = 10;
  const int nbBeads = 1000;

  std::shared_ptr<domains::Domains> allDomains(new domains::Domains(nbDomains));

  std::default_random_engine randEngine(RAND_SEED);
  std::uniform_real_distribution<double> randX(0, boxSize[0]);
  std::uniform_real_distribution<double> randY(0, boxSize[1]);
  std::uniform_real_distribution<double> randZ(0, boxSize[2]);

  cutoffgrid::CoGrid<double, cutoffgrid::CoSparseGridContainer> grid(
      coRadius, boxSize, nbDomains);

  for (int idxDom = 0; idxDom < nbDomains; ++idxDom) {
    util::Array<double> current_xyz(nbBeads, 3);
    for (int idxAtom = 0; idxAtom < nbBeads; ++idxAtom) {
      current_xyz(idxAtom, 0) = randX(randEngine);
      current_xyz(idxAtom, 1) = randY(randEngine);
      current_xyz(idxAtom, 2) = randZ(randEngine);
    }

    (*allDomains)[idxDom].reset(
        new TestingDomain<
            cutoffgrid::CoGrid<double, cutoffgrid::CoSparseGridContainer>>(
            nbBeads, idxDom, grid));
    (*allDomains)[idxDom]->setXyz(current_xyz);
    reinterpret_cast<TestingDomain<
        cutoffgrid::CoGrid<double, cutoffgrid::CoSparseGridContainer>> *>(
        (*allDomains)[idxDom].get())
        ->setCoMode(false);
  }

  const energy::ForceField forcefield = dummy_forcefield(1, -1, 1, 1, 1);

  pairkernels::PairKernelBuilder::RegisterExtraBuild(
      "TestingDomain",
      TestingPairKernel<cutoffgrid::CoGrid<
          double, cutoffgrid::CoSparseGridContainer>>::BuildTestingPairKernel);
  std::vector<std::array<std::string, 3>> defaultKernels;
  defaultKernels.emplace_back(
      std::array<std::string, 3>{{"*", "*", "TestingDomain"}});
  // We use our own kernel inside the testing class
  pairkernels::PairKernelManager kernels(defaultKernels, *allDomains);

  cutoffgrid::CutoffInteractions<double, cutoffgrid::CoSparseGridContainer>
      coalgo(coRadius, boxSize, allDomains);
  FullInteractions<double> fullalgo(boxSize, allDomains);

  // Compare the compute all

  energy::rEnergyMatrix fullRes(nbDomains, nbDomains,
                                kernels.getNbContributions());
  fullalgo.computeAll(boxSize, forcefield, kernels, fullRes);
  std::vector<int> fullInteractionsCounter(nbDomains);
  for (int idxDom = 0; idxDom < nbDomains; ++idxDom) {
    fullInteractionsCounter[idxDom] = reinterpret_cast<TestingDomain<
        cutoffgrid::CoGrid<double, cutoffgrid::CoSparseGridContainer>> *>(
                                          (*allDomains)[idxDom].get())
                                          ->getNbInteractions();
    reinterpret_cast<TestingDomain<
        cutoffgrid::CoGrid<double, cutoffgrid::CoSparseGridContainer>> *>(
        (*allDomains)[idxDom].get())
        ->resetNbInteractions();
    reinterpret_cast<TestingDomain<
        cutoffgrid::CoGrid<double, cutoffgrid::CoSparseGridContainer>> *>(
        (*allDomains)[idxDom].get())
        ->setCoMode(true);
  }

  energy::rEnergyMatrix coRes(nbDomains, nbDomains,
                              kernels.getNbContributions());

  coalgo.computeAll(boxSize, forcefield, kernels, coRes);
  for (int idxDom = 0; idxDom < nbDomains; ++idxDom) {
    EXPECT_EQ(
        (reinterpret_cast<TestingDomain<
             cutoffgrid::CoGrid<double, cutoffgrid::CoSparseGridContainer>> *>(
             (*allDomains)[idxDom].get())
             ->getNbInteractions()),
        fullInteractionsCounter[idxDom]);
    reinterpret_cast<TestingDomain<
        cutoffgrid::CoGrid<double, cutoffgrid::CoSparseGridContainer>> *>(
        (*allDomains)[idxDom].get())
        ->resetNbInteractions();
    reinterpret_cast<TestingDomain<
        cutoffgrid::CoGrid<double, cutoffgrid::CoSparseGridContainer>> *>(
        (*allDomains)[idxDom].get())
        ->setCoMode(false);
  }

  for (int idxDomTgt = 0; idxDomTgt < nbDomains; ++idxDomTgt) {
    for (int idxDomSrc = 0; idxDomSrc < nbDomains; ++idxDomSrc) {
      MY_EXPECT_DOUBLE_EQ(fullRes.getEnergy(idxDomTgt, idxDomSrc),
                          coRes.getEnergy(idxDomTgt, idxDomSrc), 1E-12);
    }
  }

  // Compare the compute one

  energy::rEnergyMatrix fullCODres(1, nbDomains, kernels.getNbContributions());
  fullalgo.computeForOneDomain((*allDomains)[0]->id(), boxSize, forcefield,
                               kernels, fullCODres);
  for (int idxDom = 0; idxDom < nbDomains; ++idxDom) {
    fullInteractionsCounter[idxDom] = reinterpret_cast<TestingDomain<
        cutoffgrid::CoGrid<double, cutoffgrid::CoSparseGridContainer>> *>(
                                          (*allDomains)[idxDom].get())
                                          ->getNbInteractions();
    reinterpret_cast<TestingDomain<
        cutoffgrid::CoGrid<double, cutoffgrid::CoSparseGridContainer>> *>(
        (*allDomains)[idxDom].get())
        ->resetNbInteractions();
    reinterpret_cast<TestingDomain<
        cutoffgrid::CoGrid<double, cutoffgrid::CoSparseGridContainer>> *>(
        (*allDomains)[idxDom].get())
        ->setCoMode(true);
  }

  energy::rEnergyMatrix coCODres(1, nbDomains, kernels.getNbContributions());
  coalgo.computeForOneDomain((*allDomains)[0]->id(), boxSize, forcefield,
                             kernels, coCODres);
  for (int idxDom = 0; idxDom < nbDomains; ++idxDom) {
    EXPECT_EQ(
        (reinterpret_cast<TestingDomain<
             cutoffgrid::CoGrid<double, cutoffgrid::CoSparseGridContainer>> *>(
             (*allDomains)[idxDom].get())
             ->getNbInteractions()),
        fullInteractionsCounter[idxDom]);
    reinterpret_cast<TestingDomain<
        cutoffgrid::CoGrid<double, cutoffgrid::CoSparseGridContainer>> *>(
        (*allDomains)[idxDom].get())
        ->resetNbInteractions();
  }

  for (int idxDomTgt = 0; idxDomTgt < nbDomains; ++idxDomTgt) {
    MY_EXPECT_DOUBLE_EQ(fullCODres.getEnergy(0, idxDomTgt),
                        coCODres.getEnergy(0, idxDomTgt), 1E-12);
  }
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

TEST(CUTOFF_TEST, testInteractionAlgorithmsOneCell) {
  const double coRadius = 10;
  const auto boxSize = util::rvec(100, 200, 300);
  const int nbDomains = 10;
  const int nbBeads = 100;

  std::shared_ptr<domains::Domains> allDomains(new domains::Domains(nbDomains));

  std::default_random_engine randEngine(RAND_SEED);
  std::uniform_real_distribution<double> randX(0, coRadius);
  std::uniform_real_distribution<double> randY(0, coRadius);
  std::uniform_real_distribution<double> randZ(0, coRadius);

  cutoffgrid::CoGrid<double, cutoffgrid::CoSparseGridContainer> grid(
      coRadius, boxSize, nbDomains);

  for (int idxDom = 0; idxDom < nbDomains; ++idxDom) {
    util::Array<double> current_xyz(nbBeads, 3);
    for (int idxAtom = 0; idxAtom < nbBeads; ++idxAtom) {
      current_xyz(idxAtom, 0) = randX(randEngine);
      current_xyz(idxAtom, 1) = randY(randEngine);
      current_xyz(idxAtom, 2) = randZ(randEngine);
    }

    (*allDomains)[idxDom].reset(
        new TestingDomain<
            cutoffgrid::CoGrid<double, cutoffgrid::CoSparseGridContainer>>(
            nbBeads, idxDom, grid));
    (*allDomains)[idxDom]->setXyz(current_xyz);
    reinterpret_cast<TestingDomain<
        cutoffgrid::CoGrid<double, cutoffgrid::CoSparseGridContainer>> *>(
        (*allDomains)[idxDom].get())
        ->setCoMode(false);
  }

  const energy::ForceField forcefield = dummy_forcefield(1, -1, 1, 1, 1);

  pairkernels::PairKernelBuilder::RegisterExtraBuild(
      "TestingDomain",
      TestingPairKernel<cutoffgrid::CoGrid<
          double, cutoffgrid::CoSparseGridContainer>>::BuildTestingPairKernel);
  std::vector<std::array<std::string, 3>> defaultKernels;
  defaultKernels.emplace_back(
      std::array<std::string, 3>{{"*", "*", "TestingDomain"}});
  // We use our own kernel inside the testing class
  pairkernels::PairKernelManager kernels(defaultKernels, *allDomains);

  cutoffgrid::CutoffInteractions<double, cutoffgrid::CoSparseGridContainer>
      coalgo(coRadius, boxSize, allDomains);
  FullInteractions<double> fullalgo(boxSize, allDomains);

  // Compare the compute all

  energy::rEnergyMatrix fullRes(nbDomains, nbDomains,
                                kernels.getNbContributions());
  fullalgo.computeAll(boxSize, forcefield, kernels, fullRes);
  for (int idxDom = 0; idxDom < nbDomains; ++idxDom) {
    const int nbInteractionsForDom = reinterpret_cast<TestingDomain<
        cutoffgrid::CoGrid<double, cutoffgrid::CoSparseGridContainer>> *>(
                                         (*allDomains)[idxDom].get())
                                         ->getNbInteractions();
    EXPECT_EQ(nbBeads * (nbBeads * (nbDomains - 1)), nbInteractionsForDom);
    reinterpret_cast<TestingDomain<
        cutoffgrid::CoGrid<double, cutoffgrid::CoSparseGridContainer>> *>(
        (*allDomains)[idxDom].get())
        ->resetNbInteractions();
    reinterpret_cast<TestingDomain<
        cutoffgrid::CoGrid<double, cutoffgrid::CoSparseGridContainer>> *>(
        (*allDomains)[idxDom].get())
        ->setCoMode(true);
  }

  energy::rEnergyMatrix coRes(nbDomains, nbDomains,
                              kernels.getNbContributions());

  coalgo.computeAll(boxSize, forcefield, kernels, coRes);
  for (int idxDom = 0; idxDom < nbDomains; ++idxDom) {
    const int nbInteractionsForDom = reinterpret_cast<TestingDomain<
        cutoffgrid::CoGrid<double, cutoffgrid::CoSparseGridContainer>> *>(
                                         (*allDomains)[idxDom].get())
                                         ->getNbInteractions();
    EXPECT_EQ(nbBeads * (nbBeads * (nbDomains - 1)), nbInteractionsForDom);
    reinterpret_cast<TestingDomain<
        cutoffgrid::CoGrid<double, cutoffgrid::CoSparseGridContainer>> *>(
        (*allDomains)[idxDom].get())
        ->resetNbInteractions();
    reinterpret_cast<TestingDomain<
        cutoffgrid::CoGrid<double, cutoffgrid::CoSparseGridContainer>> *>(
        (*allDomains)[idxDom].get())
        ->setCoMode(false);
  }

  for (int idxDomTgt = 0; idxDomTgt < nbDomains; ++idxDomTgt) {
    for (int idxDomSrc = 0; idxDomSrc < nbDomains; ++idxDomSrc) {
      // Must be equal here because no interactions are cut
      MY_EXPECT_DOUBLE_EQ(fullRes.getEnergy(idxDomTgt, idxDomSrc),
                          coRes.getEnergy(idxDomTgt, idxDomSrc), 1E-12);
    }
  }

  // Compare the compute one

  energy::rEnergyMatrix fullCODres(1, nbDomains, kernels.getNbContributions());
  fullalgo.computeForOneDomain((*allDomains)[0]->id(), boxSize, forcefield,
                               kernels, fullCODres);
  for (int idxDom = 0; idxDom < nbDomains; ++idxDom) {
    const int nbInteractionsForDom = reinterpret_cast<TestingDomain<
        cutoffgrid::CoGrid<double, cutoffgrid::CoSparseGridContainer>> *>(
                                         (*allDomains)[idxDom].get())
                                         ->getNbInteractions();
    if (idxDom == 0) {
      EXPECT_EQ(nbBeads * (nbBeads * (nbDomains - 1)), nbInteractionsForDom);
    } else {
      EXPECT_EQ(0, nbInteractionsForDom);
    }
    reinterpret_cast<TestingDomain<
        cutoffgrid::CoGrid<double, cutoffgrid::CoSparseGridContainer>> *>(
        (*allDomains)[idxDom].get())
        ->resetNbInteractions();
    reinterpret_cast<TestingDomain<
        cutoffgrid::CoGrid<double, cutoffgrid::CoSparseGridContainer>> *>(
        (*allDomains)[idxDom].get())
        ->setCoMode(true);
  }

  energy::rEnergyMatrix coCODres(1, nbDomains, kernels.getNbContributions());
  coalgo.computeForOneDomain((*allDomains)[0]->id(), boxSize, forcefield,
                             kernels, coCODres);
  for (int idxDom = 0; idxDom < nbDomains; ++idxDom) {
    const int nbInteractionsForDom = reinterpret_cast<TestingDomain<
        cutoffgrid::CoGrid<double, cutoffgrid::CoSparseGridContainer>> *>(
                                         (*allDomains)[idxDom].get())
                                         ->getNbInteractions();
    if (idxDom == 0) {
      EXPECT_EQ(nbBeads * (nbBeads * (nbDomains - 1)), nbInteractionsForDom);
    } else {
      EXPECT_EQ(0, nbInteractionsForDom);
    }
    reinterpret_cast<TestingDomain<
        cutoffgrid::CoGrid<double, cutoffgrid::CoSparseGridContainer>> *>(
        (*allDomains)[idxDom].get())
        ->resetNbInteractions();
  }

  for (int idxDomTgt = 0; idxDomTgt < nbDomains; ++idxDomTgt) {
    // Must be strictly equal here because no interactions are cut
    MY_EXPECT_DOUBLE_EQ(fullCODres.getEnergy(0, idxDomTgt),
                        coCODres.getEnergy(0, idxDomTgt), 1E-12);
  }
}
