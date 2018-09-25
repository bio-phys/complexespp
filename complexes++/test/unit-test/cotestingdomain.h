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
#ifndef COTESTINGDOMAIN_H
#define COTESTINGDOMAIN_H

#include "domains/abstractdomain.h"
#include "gtest/gtest.h"

// This domain compute the interaction given in interval
// and count them
// it does apply a manual cutoff!
// so the interaction counter and energy are not stricly equal
// to a full interaction approach (but rather to a cutoff one)
template <class GirdClass>
class TestingDomain : public domains::AbstractDomain {
  //< It should be be put in a real class! since it is managed by the cutoff
  // algo.
  const GirdClass &m_grid;
  double m_maxDistanceBetweenBeads;
  mutable int m_nbInteractions;
  bool m_ensureAreCoInteractions;

public:
  void serialize(io::Serializer &serializer) const final {
    domains::AbstractDomain::serializeCore(serializer);
  }

  TestingDomain(const int inNbAtoms, const int inId, const GirdClass &inGrid)
      : domains::AbstractDomain(
            "TestingDomain", 0, inId,
            std::vector<domains::Bead>(inNbAtoms, domains::Bead(0)),
            std::vector<double>(inNbAtoms, 0),
            std::vector<domains::BeadChainID>(inNbAtoms,
                                              domains::BeadChainID("A", 0)),
            domains::Connections(0)),
        m_grid(inGrid),
        m_maxDistanceBetweenBeads(
            2 *
            sqrt(m_grid.getRadiusPerDim()[0] * m_grid.getRadiusPerDim()[0] +
                 m_grid.getRadiusPerDim()[1] * m_grid.getRadiusPerDim()[1] +
                 m_grid.getRadiusPerDim()[2] * m_grid.getRadiusPerDim()[2])),
        m_nbInteractions(0), m_ensureAreCoInteractions(false) {}

  virtual std::string type() const override { return "TestingDomain"; }

  int getNbInteractions() const { return m_nbInteractions; }

  void setCoMode(const bool inMode) { m_ensureAreCoInteractions = inMode; }

  void resetNbInteractions() { m_nbInteractions = 0; }

  template <class ViewType>
  void energyForIntervals(ViewType inView, const domains::AbstractDomain &other,
                          const util::rvec &box,
                          const energy::ForceField &forcefield,
                          const std::pair<int, int> &myInterval,
                          const std::pair<int, int> &otherInterval,
                          const int offsetContributions) const {
    UNUSED(forcefield);
    UNUSED(inView);
    UNUSED(offsetContributions);
    if (id() == other.id()) {
      return;
    }

    const auto &xyzOther = other.xyz();
    const auto &xyzThis = xyz();
    auto d = util::rvec();

    for (auto i = myInterval.first; i < myInterval.second; ++i) {
      for (auto j = otherInterval.first; j < otherInterval.second; ++j) {
        const auto myBox = m_grid.getCoordFromPosition(
            xyzThis(i, 0), xyzThis(i, 1), xyzThis(i, 2));
        const auto otherBox = m_grid.getCoordFromPosition(
            xyzOther(j, 0), xyzOther(j, 1), xyzOther(j, 2));

        const bool areNeighbors =
            (std::abs(myBox[0] - otherBox[0]) <= 1 ||
             std::abs(myBox[0] - otherBox[0]) == m_grid.getGridSize()[0] - 1) &&
            (std::abs(myBox[1] - otherBox[1]) <= 1 ||
             std::abs(myBox[1] - otherBox[1]) == m_grid.getGridSize()[1] - 1) &&
            (std::abs(myBox[2] - otherBox[2]) <= 1 ||
             std::abs(myBox[2] - otherBox[2]) == m_grid.getGridSize()[2] - 1);

        // From the test we must be in NOT co interaction or the beads must be
        // neighbors
        // or the test is not valid
        EXPECT_TRUE(m_ensureAreCoInteractions == false || areNeighbors == true);

        // We compute the interactions only if beads are neighbors or if we
        // compute for all
        if (areNeighbors == true || m_ensureAreCoInteractions == true) {
          d[0] = xyzThis(i, 0) - xyzOther(j, 0);
          d[1] = xyzThis(i, 1) - xyzOther(j, 1);
          d[2] = xyzThis(i, 2) - xyzOther(j, 2);

          const auto r2 = util::pbc::DistSquare(d, box);
          EXPECT_LE(r2, m_maxDistanceBetweenBeads * m_maxDistanceBetweenBeads);

#pragma omp atomic update
          m_nbInteractions += 1;
        }
      }
    }
  }

  double energyForConnections(const AbstractDomain &other,
                              const util::rvec &box,
                              const energy::ForceField &forcefield) const {
    UNUSED(other);
    UNUSED(box);
    UNUSED(forcefield);
    return 0;
  }

  void energyForAllConnections(const domains::Domains &other,
                               const util::rvec &box,
                               const energy::ForceField &forcefield,
                               std::vector<double> &res) const {
    UNUSED(other);
    UNUSED(box);
    UNUSED(forcefield);
    UNUSED(res);
  }

  domains::MovedDomain move(const util::rvec &box,
                            util::RNGEngine &rng) const final {
    UNUSED(box);
    UNUSED(rng);
    return domains::MovedDomain::Fail();
  }

  std::unique_ptr<domains::AbstractDomain> copy() const final {
    return std::unique_ptr<domains::AbstractDomain>();
  }
};

#include "pairkernels/abstractpairkernel.h"

template <class GirdClass>
class TestingPairKernel : public pairkernels::AbstractPairKernel {
public:
  static std::unique_ptr<TestingPairKernel> BuildTestingPairKernel() {
    return std::make_unique<TestingPairKernel>();
  }

  explicit TestingPairKernel() : pairkernels::AbstractPairKernel(1) {}

  void compute(energy::EnergyMatrixBuffer<>::Accesser inView,
               const domains::AbstractDomain &inDom1,
               const domains::AbstractDomain &inDom2, const util::rvec &box,
               const energy::ForceField &forceField,
               const std::pair<int, int> &inInterval1,
               const std::pair<int, int> &inInterval2) const final {
    reinterpret_cast<const TestingDomain<GirdClass> *>(&inDom1)
        ->energyForIntervals(inView, inDom2, box, forceField, inInterval1,
                             inInterval2, getOffsetContributions());
  }

  static std::string GetName() { return "TestingPairKernel"; }

  std::string type() const override { return GetName(); }
};

#endif
