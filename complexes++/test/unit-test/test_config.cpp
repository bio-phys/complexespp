// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#include <unistd.h>

#include "complexes_test.h"
#include "setup/config.h"
#include "gtest/gtest.h"

const std::string configFile = dataDir + "conf.yaml";

TEST(CONFIG_TEST, read_config) {
  auto conf = setup::Config(configFile);
  EXPECT_EQ(300, conf.value<int>("montecarlo.algorithm-params.temperatur"));
}

TEST(CONFIG_TEST, read_bad_type) {
  auto conf = setup::Config(configFile);
  EXPECT_THROW_MESSAGE(conf.value<int>("montecarlo.seed"),
                       std::invalid_argument,
                       "key (montecarlo.seed) is not of type: int");
}

TEST(CONFIG_TEST, read_config_non_existing_file) {
  EXPECT_THROW(setup::Config("foo-bar-bulb-1234"), std::invalid_argument);
}
