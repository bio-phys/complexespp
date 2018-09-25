// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#include "complexes_test.h"
#include "gtest/gtest.h"
#include <cmath>

#include "domains/rigiddomain.h"
#include "pairkernels/pairkernelmanager.h"

const auto SAMETYPE = true;
const auto DIFFERENTTYPE = !SAMETYPE;

domains::Domains genDomains(const std::vector<std::string> names,
                            const bool bSameType) {
  domains::Domains topology(names.size());

  const auto Dtrans = util::rvec(.2, .2, .2);
  const auto phi = 2.0;
  for (auto idx = 0u; idx < names.size(); ++idx) {
    topology[idx] = std::make_unique<domains::Rigid>(
        names[idx], bSameType ? 0 : idx, 1,
        std::vector<domains::Bead>{domains::Bead(0)}, std::vector<double>{0},
        std::vector<domains::BeadChainID>{domains::BeadChainID("A", 1)},
        domains::Connections(0), Dtrans, phi);
  }
  return topology;
}

TEST(PAIRKERNELS, testWrongTypeId) {
  const auto topology = genDomains({{"A", "B", "C", "D", "E"}}, SAMETYPE);
  std::vector<std::array<std::string, 3>> inAllPairs;
  inAllPairs.emplace_back(std::array<std::string, 3>{{"*", "*", "None"}});

  EXPECT_THROW_MESSAGE(
      pairkernels::PairKernelManager kernels(inAllPairs, topology),
      std::invalid_argument,
      "Some domains are from the different type but have the same type ID");
}

TEST(PAIRKERNELS, testWrongType) {
  const auto topology = genDomains({{"A", "A", "A", "A", "A"}}, DIFFERENTTYPE);
  std::vector<std::array<std::string, 3>> inAllPairs;
  inAllPairs.emplace_back(std::array<std::string, 3>{{"*", "*", "None"}});

  EXPECT_THROW_MESSAGE(
      pairkernels::PairKernelManager kernels(inAllPairs, topology),
      std::invalid_argument,
      "Some domains are from the same type but have different type ID");
}

TEST(PAIRKERNELS, testAllall) {
  const auto topology = genDomains({{"A", "B", "C", "D", "E"}}, DIFFERENTTYPE);
  {
    std::vector<std::array<std::string, 3>> inAllPairs;
    inAllPairs.emplace_back(std::array<std::string, 3>{{"*", "*", "None"}});

    pairkernels::PairKernelManager kernels(inAllPairs, topology);

    for (int idx1 = 0; idx1 < 5; ++idx1) {
      for (int idx2 = 0; idx2 < 5; ++idx2) {
        EXPECT_STREQ(
            kernels
                .getKernel(topology[idx1]->typeId(), topology[idx2]->typeId())
                .type()
                .c_str(),
            "None");
      }
    }
  }
  {
    std::vector<std::array<std::string, 3>> inAllPairs;
    inAllPairs.emplace_back(
        std::array<std::string, 3>{{"default", "*", "None"}});

    pairkernels::PairKernelManager kernels(inAllPairs, topology);

    for (int idx1 = 0; idx1 < 5; ++idx1) {
      for (int idx2 = 0; idx2 < 5; ++idx2) {
        EXPECT_STREQ(
            kernels
                .getKernel(topology[idx1]->typeId(), topology[idx2]->typeId())
                .type()
                .c_str(),
            "None");
      }
    }
  }
  {
    std::vector<std::array<std::string, 3>> inAllPairs;
    inAllPairs.emplace_back(
        std::array<std::string, 3>{{"default", "default", "None"}});

    pairkernels::PairKernelManager kernels(inAllPairs, topology);

    for (int idx1 = 0; idx1 < 5; ++idx1) {
      for (int idx2 = 0; idx2 < 5; ++idx2) {
        EXPECT_STREQ(
            kernels
                .getKernel(topology[idx1]->typeId(), topology[idx2]->typeId())
                .type()
                .c_str(),
            "None");
      }
    }
  }
  {
    std::vector<std::array<std::string, 3>> inAllPairs;
    inAllPairs.emplace_back(
        std::array<std::string, 3>{{"*", "default", "None"}});

    pairkernels::PairKernelManager kernels(inAllPairs, topology);

    for (int idx1 = 0; idx1 < 5; ++idx1) {
      for (int idx2 = 0; idx2 < 5; ++idx2) {
        EXPECT_STREQ(
            kernels
                .getKernel(topology[idx1]->typeId(), topology[idx2]->typeId())
                .type()
                .c_str(),
            "None");
      }
    }
  }
}

TEST(PAIRKERNELS, testDefaultOverride) {
  const auto topology = genDomains({{"A", "B", "C", "D", "E"}}, DIFFERENTTYPE);
  {
    std::vector<std::array<std::string, 3>> inAllPairs;
    inAllPairs.emplace_back(
        std::array<std::string, 3>{{"A", "default", "LJH"}});
    inAllPairs.emplace_back(std::array<std::string, 3>{{"*", "*", "None"}});

    pairkernels::PairKernelManager kernels(inAllPairs, topology);

    for (int idx1 = 0; idx1 < 5; ++idx1) {
      for (int idx2 = 0; idx2 < 5; ++idx2) {
        EXPECT_STREQ(
            kernels
                .getKernel(topology[idx1]->typeId(), topology[idx2]->typeId())
                .type()
                .c_str(),
            "None");
      }
    }
  }
  {
    std::vector<std::array<std::string, 3>> inAllPairs;
    inAllPairs.emplace_back(std::array<std::string, 3>{{"A", "B", "LJH"}});
    inAllPairs.emplace_back(
        std::array<std::string, 3>{{"default", "default", "None"}});

    pairkernels::PairKernelManager kernels(inAllPairs, topology);

    for (int idx1 = 0; idx1 < 5; ++idx1) {
      for (int idx2 = 0; idx2 < 5; ++idx2) {
        if ((idx1 == 0 && idx2 == 1) || (idx1 == 1 && idx2 == 0)) {
          EXPECT_STREQ(
              kernels
                  .getKernel(topology[idx1]->typeId(), topology[idx2]->typeId())
                  .type()
                  .c_str(),
              "LJH");
        } else {
          EXPECT_STREQ(
              kernels
                  .getKernel(topology[idx1]->typeId(), topology[idx2]->typeId())
                  .type()
                  .c_str(),
              "None");
        }
      }
    }
  }
  {
    std::vector<std::array<std::string, 3>> inAllPairs;
    inAllPairs.emplace_back(std::array<std::string, 3>{{"A", "B", "LJH"}});
    inAllPairs.emplace_back(
        std::array<std::string, 3>{{"default", "default", "None"}});
    inAllPairs.emplace_back(
        std::array<std::string, 3>{{"C", "B", "LJHCutoff"}});

    pairkernels::PairKernelManager kernels(inAllPairs, topology);

    for (int idx1 = 0; idx1 < 5; ++idx1) {
      for (int idx2 = 0; idx2 < 5; ++idx2) {
        if ((idx1 == 0 && idx2 == 1) || (idx1 == 1 && idx2 == 0)) {
          EXPECT_STREQ(
              kernels
                  .getKernel(topology[idx1]->typeId(), topology[idx2]->typeId())
                  .type()
                  .c_str(),
              "LJH");
        } else if ((idx1 == 1 && idx2 == 2) || (idx1 == 2 && idx2 == 1)) {
          EXPECT_STREQ(
              kernels
                  .getKernel(topology[idx1]->typeId(), topology[idx2]->typeId())
                  .type()
                  .c_str(),
              "LJHCutoff");
        } else {
          EXPECT_STREQ(
              kernels
                  .getKernel(topology[idx1]->typeId(), topology[idx2]->typeId())
                  .type()
                  .c_str(),
              "None");
        }
      }
    }
  }
  {
    std::vector<std::array<std::string, 3>> inAllPairs;
    inAllPairs.emplace_back(
        std::array<std::string, 3>{{"default", "default", "None"}});
    inAllPairs.emplace_back(std::array<std::string, 3>{{"A", "B", "LJH"}});
    inAllPairs.emplace_back(
        std::array<std::string, 3>{{"B", "default", "Repulsive"}});
    inAllPairs.emplace_back(
        std::array<std::string, 3>{{"C", "*", "LJHCutoff"}});

    pairkernels::PairKernelManager kernels(inAllPairs, topology);

    for (int idx1 = 0; idx1 < 5; ++idx1) {
      for (int idx2 = 0; idx2 < 5; ++idx2) {
        if (idx1 == 2 || idx2 == 2) {
          EXPECT_STREQ(
              kernels
                  .getKernel(topology[idx1]->typeId(), topology[idx2]->typeId())
                  .type()
                  .c_str(),
              "LJHCutoff");
        } else if ((idx1 == 0 && idx2 == 1) || (idx1 == 1 && idx2 == 0)) {
          EXPECT_STREQ(
              kernels
                  .getKernel(topology[idx1]->typeId(), topology[idx2]->typeId())
                  .type()
                  .c_str(),
              "LJH");
        } else if (idx1 == 1 || idx2 == 1) {
          EXPECT_STREQ(
              kernels
                  .getKernel(topology[idx1]->typeId(), topology[idx2]->typeId())
                  .type()
                  .c_str(),
              "Repulsive");
        } else {
          EXPECT_STREQ(
              kernels
                  .getKernel(topology[idx1]->typeId(), topology[idx2]->typeId())
                  .type()
                  .c_str(),
              "None");
        }
      }
    }
  }
}

TEST(PAIRKERNELS, testFails) {
  const auto topology = genDomains({{"A", "B", "C", "D", "E"}}, DIFFERENTTYPE);
  {
    std::vector<std::array<std::string, 3>> inAllPairs;
    inAllPairs.emplace_back(std::array<std::string, 3>{{"A", "B", "None"}});

    try {
      pairkernels::PairKernelManager kernels(inAllPairs, topology);
      EXPECT_FALSE(true);  // should never be here
    } catch (std::invalid_argument& e) {
      EXPECT_TRUE(std::string(e.what()).find("PairKernelManager -- pair") !=
                      std::string::npos &&
                  std::string(e.what()).find("has never been defined") !=
                      std::string::npos);
    }
  }
  {
    std::vector<std::array<std::string, 3>> inAllPairs;
    inAllPairs.emplace_back(
        std::array<std::string, 3>{{"default", "default", "None"}});
    inAllPairs.emplace_back(std::array<std::string, 3>{{"A", "*", "None"}});
    inAllPairs.emplace_back(std::array<std::string, 3>{{"A", "B", "None"}});

    EXPECT_THROW_MESSAGE(
        pairkernels::PairKernelManager kernels(inAllPairs, topology),
        std::invalid_argument,
        "PairKernelManager -- pair A / B has already been defined");
  }
  {
    std::vector<std::array<std::string, 3>> inAllPairs;
    inAllPairs.emplace_back(
        std::array<std::string, 3>{{"*", "default", "None"}});
    inAllPairs.emplace_back(std::array<std::string, 3>{{"A", "B", "None"}});

    EXPECT_THROW_MESSAGE(
        pairkernels::PairKernelManager kernels(inAllPairs, topology),
        std::invalid_argument,
        "PairKernelManager -- pair A / B has already been defined");
  }
  {
    std::vector<std::array<std::string, 3>> inAllPairs;
    inAllPairs.emplace_back(std::array<std::string, 3>{{"A", "B", "None"}});
    inAllPairs.emplace_back(
        std::array<std::string, 3>{{"*", "default", "None"}});

    EXPECT_THROW_MESSAGE(
        pairkernels::PairKernelManager kernels(inAllPairs, topology),
        std::invalid_argument,
        "PairKernelManager -- pair A / B has already been defined");
  }
  {
    std::vector<std::array<std::string, 3>> inAllPairs;
    inAllPairs.emplace_back(std::array<std::string, 3>{{"A", "B", "None"}});
    inAllPairs.emplace_back(std::array<std::string, 3>{{"B", "A", "None"}});

    EXPECT_THROW_MESSAGE(
        pairkernels::PairKernelManager kernels(inAllPairs, topology),
        std::invalid_argument,
        "PairKernelManager -- pair B / A has already been defined");
  }
  {
    std::vector<std::array<std::string, 3>> inAllPairs;
    inAllPairs.emplace_back(std::array<std::string, 3>{
        {"default", "default", "something that does not exist"}});

    EXPECT_THROW_MESSAGE_CONTAINS(
        pairkernels::PairKernelManager kernels(inAllPairs, topology),
        std::runtime_error,
        "BuildPairKernel -- Pair kernel something that does not "
        "exist cannot be found, ");
  }
}
