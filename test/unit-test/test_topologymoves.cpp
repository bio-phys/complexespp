// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#include "gtest/gtest.h"

#include "complexes_test.h"

#include "domains/rigiddomain.h"
#include "domains/topology.h"
#include "domains/topologymoves.h"
#include "util/linalg.h"

double getCentroidScore(const domains::Domains& domains, util::rvec center) {
  double score = 0;
  for (auto& dom : domains) {
    const util::rvec domCentroid = util::centroid(dom->xyz());
    score += (center[0] - domCentroid[0]) * (center[0] - domCentroid[0]) +
             (center[1] - domCentroid[1]) * (center[1] - domCentroid[1]) +
             (center[2] - domCentroid[2]) * (center[2] - domCentroid[2]);
  }
  return score;
}

util::rvec getCentroidOfDomains(const domains::Domains& domains) {
  util::rvec centroid;
  for (auto& dom : domains) {
    centroid += util::centroid(dom->xyz());
  }
  for (int idxDim = 0; idxDim < 3; ++idxDim) {
    centroid[idxDim] /= static_cast<double>(domains.size());
  }
  return centroid;
}

domains::Domains generateDomainsWithConnections(
    const int nbDomains, const util::rvec boxSize, const int nbBeads,
    const std::vector<std::vector<int>>& connectionsMapping,
    util::RNGEngine& inRng) {
  auto dist = std::uniform_real_distribution<double>{0, 1};

  domains::Domains domains(nbDomains);

  std::vector<domains::Bead> beads(nbBeads);
  std::vector<domains::BeadChainID> beadsIds(nbBeads,
                                             domains::BeadChainID("", 0));

  for (int idx = 0; idx < nbDomains; ++idx) {
    domains::Connections connections;
    for (const int idConnection : connectionsMapping[idx]) {
      connections.emplace_back(
          std::make_unique<domains::FlatConnection>(0, idConnection, 0));
    }

    domains[idx].reset(new domains::Rigid(
        "typename", idx, idx, beads, std::vector<double>(nbBeads), beadsIds,
        connections, util::rvec(), 0));

    util::rArray pos(nbBeads, 3);

    for (int idxDim = 0; idxDim < 3; ++idxDim) {
      pos(0, idxDim) =
          dist(inRng) * boxSize[idxDim] * 0.9 + boxSize[idxDim] * 0.05;
    }
    for (int idxPart = 1; idxPart < nbBeads; ++idxPart) {
      for (int idxDim = 0; idxDim < 3; ++idxDim) {
        pos(idxPart, idxDim) = pos(idxPart - 1, idxDim) +
                               dist(inRng) * 0.05 * boxSize[idxDim] - 0.025;
      }
    }

    domains[idx]->setXyz(pos);
  }

  return domains;
}

domains::Domains generateDomains(const int nbDomains, const util::rvec boxSize,
                                 const int nbBeads, util::RNGEngine& inRng) {
  return generateDomainsWithConnections(
      nbDomains, boxSize, nbBeads, std::vector<std::vector<int>>(nbDomains),
      inRng);
}

domains::Domains cloneAllDomains(const domains::Domains& inDomains) {
  domains::Domains copy;
  copy.reserve(inDomains.size());
  for (int idxDom = 0; idxDom < static_cast<int>(inDomains.size()); ++idxDom) {
    copy.emplace_back(inDomains[idxDom]->copy());
  }
  return copy;
}

template <class FuncType>
void testTwoDomains(const domains::Domains& inDomains1,
                    const domains::Domains& inDomains2, FuncType func) {
  ASSERT_EQ(inDomains1.size(), inDomains2.size());

  for (int idxDom = 0; idxDom < static_cast<int>(inDomains1.size()); ++idxDom) {
    const auto& xyz1 = inDomains1[idxDom]->xyz();
    const auto& xyz2 = inDomains2[idxDom]->xyz();
    ASSERT_EQ(xyz1.rows(), xyz2.rows());

    for (int idxRow = 0; idxRow < xyz1.rows(); ++idxRow) {
      util::rvec v1;
      util::rvec v2;
      for (int idxDim = 0; idxDim < 3; ++idxDim) {
        v1[idxDim] = xyz1(idxRow, idxDim);
        v2[idxDim] = xyz2(idxRow, idxDim);
      }
      func(v1, v2);
    }
  }
}

TEST(TOPOLOGYMOVES_TEST, testCopy) {
  auto rng = util::RNGEngine(RAND_SEED);
  const int nbDomains = 10;
  const int nbBeads = 10;
  domains::Domains domains =
      generateDomains(nbDomains, util::rvec{1, 1, 1}, nbBeads, rng);
  {
    domains::Domains domainsCopy = cloneAllDomains(domains);

    testTwoDomains(domains, domainsCopy,
                   [](const util::rvec& v1, const util::rvec& v2) {
                     for (int idxDim = 0; idxDim < 3; ++idxDim) {
                       EXPECT_EQ(v1[idxDim], v2[idxDim]);
                     }
                   });
  }
}

TEST(TOPOLOGYMOVES_TEST, testTranslation) {
  auto rng = util::RNGEngine(RAND_SEED);
  const int nbDomains = 10;
  const int nbBeads = 10;
  domains::Domains domains =
      generateDomains(nbDomains, util::rvec{1, 1, 1}, nbBeads, rng);
  std::vector<int> domainsIds = util::arange<int>(nbDomains);
  {
    const double translationCoef = 1;
    domains::Domains domainsCopy = cloneAllDomains(domains);
    domains::topologyMoves::Translation(domainsIds, translationCoef,
                                        &domainsCopy, &rng);

    testTwoDomains(
        domains, domainsCopy,
        [translationCoef](const util::rvec& v1, const util::rvec& v2) {
          for (int idxDim = 0; idxDim < 3; ++idxDim) {
            EXPECT_TRUE(v1[idxDim] - translationCoef <= v2[idxDim]);
            EXPECT_TRUE(v2[idxDim] <= v1[idxDim] + translationCoef);
          }
        });
  }
  {
    const double translationCoef = 10;
    domains::Domains domainsCopy = cloneAllDomains(domains);
    domains::topologyMoves::Translation(domainsIds, translationCoef,
                                        &domainsCopy, &rng);

    testTwoDomains(
        domains, domainsCopy,
        [translationCoef](const util::rvec& v1, const util::rvec& v2) {
          for (int idxDim = 0; idxDim < 3; ++idxDim) {
            EXPECT_TRUE(v1[idxDim] - translationCoef <= v2[idxDim]);
            EXPECT_TRUE(v2[idxDim] <= v1[idxDim] + translationCoef);
          }
        });
  }
  {
    const double translationCoef = 0;
    domains::Domains domainsCopy = cloneAllDomains(domains);
    domains::topologyMoves::Translation(domainsIds, translationCoef,
                                        &domainsCopy, &rng);

    testTwoDomains(domains, domainsCopy,
                   [](const util::rvec& v1, const util::rvec& v2) {
                     for (int idxDim = 0; idxDim < 3; ++idxDim) {
                       EXPECT_EQ(v1[idxDim], v2[idxDim]);
                     }
                   });
  }
}

TEST(TOPOLOGYMOVES_TEST, testCompactSame) {
  auto rng = util::RNGEngine(RAND_SEED);
  const int nbDomains = 10;
  const int nbBeads = 10;
  const util::rvec boxSize(1, 1, 1);

  std::vector<std::vector<int>> connectionIds(nbDomains);
  for (int idxDom = 0; idxDom < nbDomains; ++idxDom) {
    connectionIds[idxDom].push_back((idxDom + 1) % nbDomains);
  }

  domains::Domains domains = generateDomainsWithConnections(
      nbDomains, boxSize, nbBeads, connectionIds, rng);
  std::vector<int> domainsIds = util::arange<int>(nbDomains);
  {
    domains::Domains domainsCopy = cloneAllDomains(domains);

    const auto com =
        domains::topologyMoves::CompactTopology(boxSize, domainsIds, &domains);
    const util::rvec computedCentroid = getCentroidOfDomains(domains);
    for (int idxDim = 0; idxDim < 3; ++idxDim) {
      MY_EXPECT_DOUBLE_EQ(com[idxDim], computedCentroid[idxDim], 1E-13);
    }

    const auto comCopy = domains::topologyMoves::CompactTopology(
        boxSize, domainsIds, &domainsCopy);
    const util::rvec computedCentroidCopy = getCentroidOfDomains(domainsCopy);
    for (int idxDim = 0; idxDim < 3; ++idxDim) {
      MY_EXPECT_DOUBLE_EQ(comCopy[idxDim], computedCentroidCopy[idxDim], 1E-13);
    }

    for (int idxDim = 0; idxDim < 3; ++idxDim) {
      MY_EXPECT_DOUBLE_EQ(com[idxDim], comCopy[idxDim], 1E-13);
    }

    testTwoDomains(domains, domainsCopy,
                   [](const util::rvec& v1, const util::rvec& v2) {
                     for (int idxDim = 0; idxDim < 3; ++idxDim) {
                       MY_EXPECT_DOUBLE_EQ(v1[idxDim], v2[idxDim], 1E-13);
                     }
                   });
  }
}

TEST(TOPOLOGYMOVES_TEST, testCompactPbc) {
  auto rng = util::RNGEngine(RAND_SEED);
  const int nbDomains = 10;
  const int nbBeads = 10;
  const util::rvec boxSize(1, 1, 1);

  std::vector<std::vector<int>> connectionIds(nbDomains);
  for (int idxDom = 0; idxDom < nbDomains; ++idxDom) {
    connectionIds[idxDom].push_back((idxDom + 1) % nbDomains);
  }

  domains::Domains domains = generateDomainsWithConnections(
      nbDomains, boxSize, nbBeads, connectionIds, rng);
  std::vector<int> domainsIds = util::arange<int>(nbDomains);

  {
    domains::Domains domainsCopy = cloneAllDomains(domains);
    auto comCopy = domains::topologyMoves::CompactTopology(boxSize, domainsIds,
                                                           &domainsCopy);

    const util::rvec computedCentroid = getCentroidOfDomains(domainsCopy);
    for (int idxDim = 0; idxDim < 3; ++idxDim) {
      MY_EXPECT_DOUBLE_EQ(comCopy[idxDim], computedCentroid[idxDim], 1E-13);
    }

    testTwoDomains(domains, domainsCopy, [&boxSize](const util::rvec& v1,
                                                    const util::rvec& v2) {
      for (int idxDim = 0; idxDim < 3; ++idxDim) {
        const double v1tobox = util::pbc::ToBox(v1[idxDim], boxSize[idxDim]);
        const double v2tobox = util::pbc::ToBox(v2[idxDim], boxSize[idxDim]);
        MY_EXPECT_DOUBLE_EQ(v1tobox, v2tobox, 1E-11);
      }
    });
  }
}

TEST(TOPOLOGYMOVES_TEST, testCompactSameQuality) {
  auto rng = util::RNGEngine(RAND_SEED);
  const int nbDomains = 10;
  const int nbBeads = 10;
  const util::rvec boxSize(1, 1, 1);

  std::vector<std::vector<int>> connectionIds(nbDomains);
  for (int idxDom = 0; idxDom < nbDomains; ++idxDom) {
    connectionIds[idxDom].push_back((idxDom + 1) % nbDomains);
  }

  domains::Domains domains = generateDomainsWithConnections(
      nbDomains, boxSize, nbBeads, connectionIds, rng);

  std::vector<int> domainsIds = util::arange<int>(nbDomains);

  {
    domains::Domains domainsCopy = cloneAllDomains(domains);

    for (auto& dom : domains) {
      util::translateByConstantInPlace(boxSize[0] * 10, &dom->xyz());
    }

    const auto com =
        domains::topologyMoves::CompactTopology(boxSize, domainsIds, &domains);
    const auto comCopy = domains::topologyMoves::CompactTopology(
        boxSize, domainsIds, &domainsCopy);

    const double score = getCentroidScore(domains, com);
    const double scoreCopy = getCentroidScore(domainsCopy, comCopy);
    MY_EXPECT_DOUBLE_EQ(score, scoreCopy, 1E-13);

    for (int idxDim = 0; idxDim < 3; ++idxDim) {
      const double com1tobox = util::pbc::ToBox(com[idxDim], boxSize[idxDim]);
      const double com2tobox =
          util::pbc::ToBox(comCopy[idxDim], boxSize[idxDim]);
      MY_EXPECT_DOUBLE_EQ(com1tobox, com2tobox, 1E-11);
    }

    testTwoDomains(domains, domainsCopy, [&boxSize](const util::rvec& v1,
                                                    const util::rvec& v2) {
      for (int idxDim = 0; idxDim < 3; ++idxDim) {
        const double v1tobox = util::pbc::ToBox(v1[idxDim], boxSize[idxDim]);
        const double v2tobox = util::pbc::ToBox(v2[idxDim], boxSize[idxDim]);
        MY_EXPECT_DOUBLE_EQ(v1tobox, v2tobox, 1E-11);
      }
    });
  }
}

TEST(TOPOLOGYMOVES_TEST, testCompactKnown) {
  auto rng = util::RNGEngine(RAND_SEED);
  const int nbDomains = 10;
  const int nbBeads = 1;
  const util::rvec subBoxSize(0.1, 0.1, 0.1);

  std::vector<std::vector<int>> connectionIds(nbDomains);
  for (int idxDom = 0; idxDom < nbDomains; ++idxDom) {
    connectionIds[idxDom].push_back((idxDom + 1) % nbDomains);
  }

  domains::Domains domains = generateDomainsWithConnections(
      nbDomains, subBoxSize, nbBeads, connectionIds, rng);
  std::vector<int> domainsIds = util::arange<int>(nbDomains);

  const util::rvec boxSize(1, 1, 1);
  domains::Domains domainsCopy = cloneAllDomains(domains);
  for (auto& dom : domainsCopy) {
    util::translateByConstantInPlace(boxSize[0] * (dom->id() + 1), &dom->xyz());
  }

  testTwoDomains(
      domains, domainsCopy,
      [&subBoxSize, &boxSize](const util::rvec& v1, const util::rvec& v2) {
        for (int idxDim = 0; idxDim < 3; ++idxDim) {
          EXPECT_TRUE(v1[idxDim] < subBoxSize[idxDim]);
          const double v2tobox = util::pbc::ToBox(v2[idxDim], boxSize[idxDim]);
          MY_EXPECT_DOUBLE_EQ(v1[idxDim], v2tobox, 1E-11);
        }
      });

  const auto com =
      domains::topologyMoves::CompactTopology(boxSize, domainsIds, &domains);
  const auto comCopy = domains::topologyMoves::CompactTopology(
      boxSize, domainsIds, &domainsCopy);

  const double score = getCentroidScore(domains, com);
  const double scoreCopy = getCentroidScore(domainsCopy, comCopy);
  MY_EXPECT_DOUBLE_EQ(score, scoreCopy, 1E-13);

  util::rvec shift;
  for (int idxDim = 0; idxDim < 3; ++idxDim) {
    EXPECT_TRUE(com[idxDim] < subBoxSize[idxDim]);
    const double com2tobox = util::pbc::ToBox(comCopy[idxDim], boxSize[idxDim]);
    MY_EXPECT_DOUBLE_EQ(com[idxDim], com2tobox, 1E-11);
    shift[idxDim] = com[idxDim] - comCopy[idxDim];
  }

  testTwoDomains(
      domains, domainsCopy,
      [&shift, &subBoxSize](const util::rvec& v1, const util::rvec& v2) {
        for (int idxDim = 0; idxDim < 3; ++idxDim) {
          EXPECT_TRUE(v1[idxDim] < subBoxSize[idxDim]);
          EXPECT_TRUE(v2[idxDim] + shift[idxDim] < subBoxSize[idxDim]);
          MY_EXPECT_DOUBLE_EQ(v1[idxDim], v2[idxDim] + shift[idxDim], 1E-11);
        }
      });
}

TEST(TOPOLOGYMOVES_TEST, testRotation0) {
  auto rng = util::RNGEngine(RAND_SEED);
  const int nbDomains = 10;
  const int nbBeads = 10;
  const util::rvec subBoxSize(0.1, 0.1, 0.1);
  const util::rvec boxSize(1, 1, 1);

  std::vector<std::vector<int>> connectionIds(nbDomains);
  for (int idxDom = 0; idxDom < nbDomains; ++idxDom) {
    connectionIds[idxDom].push_back((idxDom + 1) % nbDomains);
  }

  domains::Domains domains = generateDomainsWithConnections(
      nbDomains, subBoxSize, nbBeads, connectionIds, rng);
  std::vector<int> domainsIds = util::arange<int>(nbDomains);
  const util::rvec computedCentroid = getCentroidOfDomains(domains);

  {
    domains::Domains domainsCopy = cloneAllDomains(domains);
    const auto comCopy = domains::topologyMoves::CompactTopology(
        boxSize, domainsIds, &domainsCopy);

    for (int idxDim = 0; idxDim < 3; ++idxDim) {
      MY_EXPECT_DOUBLE_EQ(computedCentroid[idxDim], comCopy[idxDim], 1E-13);
    }

    testTwoDomains(domains, domainsCopy,
                   [](const util::rvec& v1, const util::rvec& v2) {
                     for (int idxDim = 0; idxDim < 3; ++idxDim) {
                       EXPECT_EQ(v1[idxDim], v2[idxDim]);
                     }
                   });
  }
  {
    const double rotationCoef = 0;
    domains::Domains domainsCopy = cloneAllDomains(domains);
    domains::topologyMoves::Rotation(boxSize, domainsIds, rotationCoef,
                                     &domainsCopy, &rng);

    testTwoDomains(
        domains, domainsCopy,
        [&computedCentroid](const util::rvec& v1, const util::rvec& v2) {
          double dist1 = 0;
          double dist2 = 0;
          for (int idxDim = 0; idxDim < 3; ++idxDim) {
            dist1 += (v1[idxDim] - computedCentroid[idxDim]) *
                     (v1[idxDim] - computedCentroid[idxDim]);
            dist2 += (v2[idxDim] - computedCentroid[idxDim]) *
                     (v2[idxDim] - computedCentroid[idxDim]);
          }
          MY_EXPECT_DOUBLE_EQ(dist1, dist2, 1E-13);
        });

    testTwoDomains(domains, domainsCopy,
                   [](const util::rvec& v1, const util::rvec& v2) {
                     for (int idxDim = 0; idxDim < 3; ++idxDim) {
                       EXPECT_EQ(v1[idxDim], v2[idxDim]);
                     }
                   });
  }
  {
    const double rotationCoef = std::numeric_limits<double>::epsilon();
    domains::Domains domainsCopy = cloneAllDomains(domains);
    domains::topologyMoves::Rotation(boxSize, domainsIds, rotationCoef,
                                     &domainsCopy, &rng);

    testTwoDomains(
        domains, domainsCopy,
        [&computedCentroid](const util::rvec& v1, const util::rvec& v2) {
          double dist1 = 0;
          double dist2 = 0;
          for (int idxDim = 0; idxDim < 3; ++idxDim) {
            dist1 += (v1[idxDim] - computedCentroid[idxDim]) *
                     (v1[idxDim] - computedCentroid[idxDim]);
            dist2 += (v2[idxDim] - computedCentroid[idxDim]) *
                     (v2[idxDim] - computedCentroid[idxDim]);
          }
          MY_EXPECT_DOUBLE_EQ(dist1, dist2, 1E-13);
        });

    testTwoDomains(domains, domainsCopy,
                   [](const util::rvec& v1, const util::rvec& v2) {
                     for (int idxDim = 0; idxDim < 3; ++idxDim) {
                       MY_EXPECT_DOUBLE_EQ(v1[idxDim], v2[idxDim], 1E-13);
                     }
                   });
  }
}

TEST(TOPOLOGYMOVES_TEST, testRotation) {
  auto rng = util::RNGEngine(RAND_SEED);
  const int nbDomains = 10;
  const int nbBeads = 10;
  const util::rvec subBoxSize(0.1, 0.1, 0.1);
  const util::rvec boxSize(1, 1, 1);

  std::vector<std::vector<int>> connectionIds(nbDomains);
  for (int idxDom = 0; idxDom < nbDomains; ++idxDom) {
    connectionIds[idxDom].push_back((idxDom + 1) % nbDomains);
  }

  domains::Domains domains = generateDomainsWithConnections(
      nbDomains, subBoxSize, nbBeads, connectionIds, rng);

  std::vector<int> domainsIds = util::arange<int>(nbDomains);
  const util::rvec computedCentroid = getCentroidOfDomains(domains);

  {
    domains::Domains domainsCopy = cloneAllDomains(domains);
    const auto comCopy = domains::topologyMoves::CompactTopology(
        boxSize, domainsIds, &domainsCopy);

    for (int idxDim = 0; idxDim < 3; ++idxDim) {
      MY_EXPECT_DOUBLE_EQ(computedCentroid[idxDim], comCopy[idxDim], 1E-13);
    }

    testTwoDomains(domains, domainsCopy,
                   [](const util::rvec& v1, const util::rvec& v2) {
                     for (int idxDim = 0; idxDim < 3; ++idxDim) {
                       EXPECT_EQ(v1[idxDim], v2[idxDim]);
                     }
                   });
  }
  {
    const double rotationCoef = 1;
    domains::Domains domainsCopy = cloneAllDomains(domains);
    domains::topologyMoves::Rotation(boxSize, domainsIds, rotationCoef,
                                     &domainsCopy, &rng);

    testTwoDomains(
        domains, domainsCopy,
        [&computedCentroid](const util::rvec& v1, const util::rvec& v2) {
          double dist1 = 0;
          double dist2 = 0;
          for (int idxDim = 0; idxDim < 3; ++idxDim) {
            dist1 += (v1[idxDim] - computedCentroid[idxDim]) *
                     (v1[idxDim] - computedCentroid[idxDim]);
            dist2 += (v2[idxDim] - computedCentroid[idxDim]) *
                     (v2[idxDim] - computedCentroid[idxDim]);
          }
          MY_EXPECT_DOUBLE_EQ(dist1, dist2, 1E-13);
        });
  }
}

TEST(TOPOLOGYMOVES_TEST, testValidityTopology) {
  auto rng = util::RNGEngine(RAND_SEED);
  const int nbDomains = 10;
  const int nbBeads = 10;
  const util::rvec boxSize(1, 1, 1);
  std::vector<int> domainsIds = util::arange<int>(nbDomains);

  bool topologyIsValid;
  int badDomainId;
  {
    domains::Domains domains =
        generateDomains(nbDomains, boxSize, nbBeads, rng);

    {
      domains::Topology topology(domainsIds, false);
      std::tie(topologyIsValid, badDomainId) =
          topology.allConnectionsAreValid(domains);
      EXPECT_TRUE(topologyIsValid);
      EXPECT_EQ(badDomainId, -1);
    }
    {
      domains::Topology topology(std::vector<int>{0}, false);
      std::tie(topologyIsValid, badDomainId) =
          topology.allConnectionsAreValid(domains);
      EXPECT_TRUE(topologyIsValid);
      EXPECT_EQ(badDomainId, -1);
    }
    {
      domains::Topology topology(std::vector<int>{-1}, false);
      std::tie(topologyIsValid, badDomainId) =
          topology.allConnectionsAreValid(domains);
      EXPECT_FALSE(topologyIsValid);
      EXPECT_EQ(badDomainId, -1);
    }
    {
      domains::Topology topology(std::vector<int>{nbDomains}, false);
      std::tie(topologyIsValid, badDomainId) =
          topology.allConnectionsAreValid(domains);
      EXPECT_FALSE(topologyIsValid);
      EXPECT_EQ(badDomainId, nbDomains);
    }
  }
  {
    std::vector<std::vector<int>> connectionIds(nbDomains);
    for (int idxDom = 0; idxDom < nbDomains; ++idxDom) {
      connectionIds[idxDom].push_back((idxDom + 1) % nbDomains);
    }

    domains::Domains domains = generateDomainsWithConnections(
        nbDomains, boxSize, nbBeads, connectionIds, rng);
    domains::Topology topology(domainsIds, false);
    std::tie(topologyIsValid, badDomainId) =
        topology.allConnectionsAreValid(domains);
    EXPECT_TRUE(topologyIsValid);
    EXPECT_EQ(badDomainId, -1);
  }
  {
    std::vector<std::vector<int>> connectionIds(nbDomains);
    for (int idxDom = 0; idxDom < nbDomains; ++idxDom) {
      connectionIds[idxDom].push_back((idxDom + 1) % nbDomains);
      connectionIds[idxDom].push_back((idxDom + 2) % nbDomains);
      connectionIds[idxDom].push_back((idxDom + 3) % nbDomains);
    }

    domains::Domains domains = generateDomainsWithConnections(
        nbDomains, boxSize, nbBeads, connectionIds, rng);
    domains::Topology topology(domainsIds, false);
    std::tie(topologyIsValid, badDomainId) =
        topology.allConnectionsAreValid(domains);
    EXPECT_TRUE(topologyIsValid);
    EXPECT_EQ(badDomainId, -1);
  }
  {
    std::vector<std::vector<int>> connectionIds(nbDomains);
    for (int idxDom = 0; idxDom < nbDomains; ++idxDom) {
      connectionIds[idxDom].push_back((idxDom + 1) % nbDomains);
    }
    connectionIds[0].push_back(-1);

    domains::Domains domains = generateDomainsWithConnections(
        nbDomains, boxSize, nbBeads, connectionIds, rng);
    domains::Topology topology(domainsIds, false);
    std::tie(topologyIsValid, badDomainId) =
        topology.allConnectionsAreValid(domains);
    EXPECT_FALSE(topologyIsValid);
    EXPECT_EQ(badDomainId, 0);
  }
  {
    std::vector<std::vector<int>> connectionIds(nbDomains);
    connectionIds[nbDomains - 1].push_back(nbDomains);

    domains::Domains domains = generateDomainsWithConnections(
        nbDomains, boxSize, nbBeads, connectionIds, rng);
    domains::Topology topology(domainsIds, false);
    std::tie(topologyIsValid, badDomainId) =
        topology.allConnectionsAreValid(domains);
    EXPECT_FALSE(topologyIsValid);
    EXPECT_EQ(badDomainId, nbDomains - 1);
  }
}

TEST(TOPOLOGYMOVES_TEST, testConnectivityTopology) {
  auto rng = util::RNGEngine(RAND_SEED);
  const int nbDomains = 10;
  const int nbBeads = 10;
  const util::rvec boxSize(1, 1, 1);
  std::vector<int> domainsIds = util::arange<int>(nbDomains);

  bool topologyIsConnected;
  int badDomainId;

  {
    domains::Domains domains =
        generateDomains(nbDomains, boxSize, nbBeads, rng);

    {
      domains::Topology topology(domainsIds, false);
      std::tie(topologyIsConnected, badDomainId) =
          topology.allDomainsConnected(domains);
      EXPECT_FALSE(topologyIsConnected);
      EXPECT_EQ(badDomainId, 1);
    }
    {
      domains::Topology topology(std::vector<int>{0}, false);
      std::tie(topologyIsConnected, badDomainId) =
          topology.allDomainsConnected(domains);
      EXPECT_TRUE(topologyIsConnected);
      EXPECT_EQ(badDomainId, -1);
    }
    {
      domains::Topology topology(std::vector<int>{-1}, false);
      std::tie(topologyIsConnected, badDomainId) =
          topology.allDomainsConnected(domains);
      EXPECT_FALSE(topologyIsConnected);
      EXPECT_EQ(badDomainId, -1);
    }
    {
      domains::Topology topology(std::vector<int>{nbDomains}, false);
      std::tie(topologyIsConnected, badDomainId) =
          topology.allDomainsConnected(domains);
      EXPECT_FALSE(topologyIsConnected);
      EXPECT_EQ(badDomainId, nbDomains);
    }
  }
  {
    std::vector<std::vector<int>> connectionIds(nbDomains);
    for (int idxDom = 0; idxDom < nbDomains; ++idxDom) {
      connectionIds[idxDom].push_back((idxDom + 1) % nbDomains);
    }

    domains::Domains domains = generateDomainsWithConnections(
        nbDomains, boxSize, nbBeads, connectionIds, rng);
    domains::Topology topology(domainsIds, false);
    std::tie(topologyIsConnected, badDomainId) =
        topology.allDomainsConnected(domains);
    EXPECT_TRUE(topologyIsConnected);
    EXPECT_EQ(badDomainId, -1);
  }
  {
    std::vector<std::vector<int>> connectionIds(nbDomains);
    for (int idxDom = 0; idxDom < nbDomains; ++idxDom) {
      connectionIds[idxDom].push_back((idxDom + 1) % nbDomains);
      connectionIds[idxDom].push_back((idxDom + 2) % nbDomains);
      connectionIds[idxDom].push_back((idxDom + 3) % nbDomains);
    }

    domains::Domains domains = generateDomainsWithConnections(
        nbDomains, boxSize, nbBeads, connectionIds, rng);
    domains::Topology topology(domainsIds, false);
    std::tie(topologyIsConnected, badDomainId) =
        topology.allDomainsConnected(domains);
    EXPECT_TRUE(topologyIsConnected);
    EXPECT_EQ(badDomainId, -1);
  }
  {
    std::vector<std::vector<int>> connectionIds(nbDomains);
    for (int idxDom = 0; idxDom < nbDomains; ++idxDom) {
      connectionIds[idxDom].push_back((idxDom + 1) % nbDomains);
    }
    connectionIds[0].clear();

    domains::Domains domains = generateDomainsWithConnections(
        nbDomains, boxSize, nbBeads, connectionIds, rng);
    domains::Topology topology(domainsIds, false);
    std::tie(topologyIsConnected, badDomainId) =
        topology.allDomainsConnected(domains);
    EXPECT_FALSE(topologyIsConnected);
    EXPECT_EQ(badDomainId, 1);
  }
  {
    std::vector<std::vector<int>> connectionIds(nbDomains);
    domains::Domains domains = generateDomainsWithConnections(
        nbDomains, boxSize, nbBeads, connectionIds, rng);
    domains::Topology topology(domainsIds, false);
    std::tie(topologyIsConnected, badDomainId) =
        topology.allDomainsConnected(domains);
    EXPECT_FALSE(topologyIsConnected);
    EXPECT_EQ(badDomainId, 1);
  }
  {
    std::vector<std::vector<int>> connectionIds(nbDomains);
    connectionIds[0].push_back(1);
    domains::Domains domains = generateDomainsWithConnections(
        nbDomains, boxSize, nbBeads, connectionIds, rng);
    domains::Topology topology(std::vector<int>{0}, false);
    std::tie(topologyIsConnected, badDomainId) =
        topology.allDomainsConnected(domains);
    EXPECT_FALSE(topologyIsConnected);
    EXPECT_EQ(badDomainId, 0);
  }
}
