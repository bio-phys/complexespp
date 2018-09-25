// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#include <string>

#include "complexes_test.h"
#include "gtest/gtest.h"

#include "io/yaml.h"

TEST(YAML, LOAD) {
  const auto yaml =
      io::YamlNode(YAML::Load("0:\n"
                              "  domain-a: ['r_beginning', A 1]\n"
                              "1:\n"
                              "  params: {}\n"
                              "2:\n"
                              "  params: {}\n"
                              "3:\n"
                              "  params: {}\n"));

  auto keys = yaml.keys<int>();
  std::sort(std::begin(keys), std::end(keys));
  EXPECT_EQ(keys.size(), 4);
  for (auto i = 0u; i < keys.size(); ++i) {
    EXPECT_EQ(keys[i], i);
  }
}

TEST(YAML, UNDEFINED_NODE) {
  const auto node = io::YamlNode(YAML::Load(""));
  EXPECT_THROW_MESSAGE_CONTAINS(node["foo"], std::runtime_error,
                                "Reading key /.foo");
}

TEST(YAML, ERROR_MESSAGE) {
  const auto yaml =
      io::YamlNode(YAML::Load("foo:\n"
                              "  bazbaz: 42"));

  try {
    yaml["foo"]["bar"].as<int>();
  } catch (std::runtime_error& e) {
    // Only check for substring because actual error message depends on yaml-cpp
    // version
    const std::string expectedStr("ERROR in YAML: Reading key /.foo.bar\n  ");
    EXPECT_STREQ(expectedStr.c_str(),
                 std::string(e.what()).substr(0, expectedStr.length()).c_str());
  }
}

TEST(YAML, parseVector) {
  const auto yaml =
      io::YamlNode(YAML::Load("- 1\n"
                              "- 2\n"
                              "- 3\n"
                              "- 4\n"
                              "- 5\n"));
  const auto v = io::parseVector<std::size_t>(yaml);

  EXPECT_EQ(v.size(), 5);
  for (auto i = 0u; i < v.size(); ++i) {
    EXPECT_EQ(v[i], i + 1);
  }
}

TEST(YAML, parse3DVector) {
  const auto yaml =
      io::YamlNode(YAML::Load("- 1\n"
                              "- 2\n"
                              "- 3\n"));
  const auto v = io::parse3DVector<std::size_t>(yaml);

  for (auto i = 0; i < v.size(); ++i) {
    EXPECT_EQ(v[i], i + 1);
  }
}

TEST(YAML, parse3DVector_EXCEPTION) {
  const auto yaml =
      io::YamlNode(YAML::Load("- 1\n"
                              "- 2\n"
                              "- 3\n"
                              "- 4\n"));
  EXPECT_THROW(io::parse3DVector<int>(yaml), std::runtime_error);
}

TEST(YAML, parseArray) {
  const auto yaml =
      io::YamlNode(YAML::Load("- [1, 1, 1]\n"
                              "- [2, 2, 2]\n"
                              "- [3, 3, 3]\n"
                              "- [4, 4, 4]\n"
                              "- [5, 5, 5]\n"));
  const auto v = io::parseArray<int>(yaml);

  ASSERT_EQ(v.rows(), 5);
  ASSERT_EQ(v.cols(), 3);
  for (auto i = 0; i < v.rows(); ++i) {
    for (auto j = 0; j < v.cols(); ++j) {
      EXPECT_EQ(v(i, j), i + 1);
    }
  }
}

TEST(YAML, parseArray_EXCEPTION) {
  const auto yaml =
      io::YamlNode(YAML::Load("- [1, 1, 1]\n"
                              "- [2, 2, 2]\n"
                              "- [3, 3]\n"
                              "- [4, 4, 4]\n"
                              "- [5, 5, 5]\n"));
  EXPECT_THROW(io::parseArray<int>(yaml), std::runtime_error);
}

TEST(YAML, parseFixedArray) {
  const auto yaml =
      io::YamlNode(YAML::Load("- [1, 1, 1]\n"
                              "- [2, 2, 2]\n"
                              "- [3, 3, 3]\n"
                              "- [4, 4, 4]\n"
                              "- [5, 5, 5]\n"));
  const auto v = io::parseFixedArray<int, 5, 3>(yaml);
  for (auto i = 0; i < 5; ++i) {
    for (auto j = 0; j < 3; ++j) {
      EXPECT_EQ(v(i, j), i + 1);
    }
  }
}

TEST(YAML, parseFixedArray_EXCEPTION) {
  const auto yaml =
      io::YamlNode(YAML::Load("- [1, 1, 1]\n"
                              "- [2, 2, 2]\n"
                              "- [3, 3]\n"
                              "- [4, 4, 4]\n"
                              "- [5, 5, 5]\n"));
  auto wrong_shape = [](const io::YamlNode& node) {
    io::parseFixedArray<int, 5, 3>(node);
  };
  EXPECT_THROW(wrong_shape(yaml), std::runtime_error);

  auto wrong_cols = [](const io::YamlNode& node) {
    io::parseFixedArray<int, 5, 5>(node);
  };
  EXPECT_THROW(wrong_cols(yaml), std::runtime_error);

  auto wrong_rows = [](const io::YamlNode& node) {
    io::parseFixedArray<int, 2, 3>(node);
  };
  EXPECT_THROW(wrong_rows(yaml), std::runtime_error);
}
