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
