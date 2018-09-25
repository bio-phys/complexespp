// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#include "gtest/gtest.h"

#include "util/string.h"

TEST(STRING_TEST, splitStrTest_single_token) {
  auto const tokens = util::splitStr("foo.bar.baz", ".");
  EXPECT_EQ(3u, tokens.size());
  EXPECT_STREQ("foo", tokens[0].c_str());
  EXPECT_STREQ("bar", tokens[1].c_str());
  EXPECT_STREQ("baz", tokens[2].c_str());
}

TEST(STRING_TEST, splitStrTest_double_token) {
  auto const tokens = util::splitStr("foo  bar", " ");
  EXPECT_EQ(2u, tokens.size());
  EXPECT_STREQ("foo", tokens[0].c_str());
  EXPECT_STREQ("bar", tokens[1].c_str());
}

TEST(STRING_TEST, splitStrTest_trailing_delimiter) {
  auto const tokens = util::splitStr("foo bar  ", " ");
  EXPECT_EQ(2u, tokens.size());
  EXPECT_STREQ("foo", tokens[0].c_str());
  EXPECT_STREQ("bar", tokens[1].c_str());
}

TEST(STRING_TEST, splitStrTest_prepended_delimiter) {
  auto const tokens = util::splitStr("  foo bar", " ");
  EXPECT_EQ(2u, tokens.size());
  EXPECT_STREQ("foo", tokens[0].c_str());
  EXPECT_STREQ("bar", tokens[1].c_str());
}

TEST(STRING_TEST, truncateLeft) {
  EXPECT_EQ(1, util::truncateLeft(100001));
  EXPECT_EQ(1, util::truncateLeft(1000001));
  EXPECT_EQ(1, util::truncateLeft(1));
  EXPECT_EQ(42, util::truncateLeft(42));
  EXPECT_EQ(99999, util::truncateLeft(99999));
  EXPECT_EQ(99998, util::truncateLeft(99998));
  EXPECT_EQ(2, util::truncateLeft(102, 2));
  EXPECT_EQ(99, util::truncateLeft(99, 2));
}

int truncateLeft(const int val, const int base = 10, const int max = 99999);
