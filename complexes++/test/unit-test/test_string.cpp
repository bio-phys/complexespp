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

#include "util/string.h"

TEST(STRING_TEST, splitStrTest_single_token) {
  auto const tokens = util::splitStr("foo.bar.baz", '.');
  EXPECT_EQ(3u, tokens.size());
  EXPECT_STREQ("foo", tokens[0].c_str());
  EXPECT_STREQ("bar", tokens[1].c_str());
  EXPECT_STREQ("baz", tokens[2].c_str());
}

TEST(STRING_TEST, splitStrTest_double_token) {
  auto const tokens = util::splitStr("foo  bar", ' ');
  EXPECT_EQ(2u, tokens.size());
  EXPECT_STREQ("foo", tokens[0].c_str());
  EXPECT_STREQ("bar", tokens[1].c_str());
}

TEST(STRING_TEST, splitStrTest_trailing_delimiter) {
  auto const tokens = util::splitStr("foo bar  ", ' ');
  EXPECT_STREQ("foo", tokens[0].c_str());
  EXPECT_STREQ("bar", tokens[1].c_str());
  EXPECT_EQ(2u, tokens.size());
}

TEST(STRING_TEST, splitStrTest_prepended_delimiter) {
  auto const tokens = util::splitStr("  foo bar", ' ');
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
