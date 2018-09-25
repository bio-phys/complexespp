// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#include "complexes_test.h"
#include "gtest/gtest.h"

#include "util/array.h"
#include "util/linalg.h"
#include "util/pbc.h"
#include "util/random.h"

template <typename T>
util::vec<T> closestImageExhaustive(const util::vec<T>& a,
                                    const util::vec<T>& b,
                                    const util::vec<T>& box) {
  auto closestImage = util::vec<T>(0, 0, 0);
  const auto values = std::array<T, 3>{{-1, 0, 1}};
  auto shortest = 1e10;
  for (const auto i : values) {
    for (const auto j : values) {
      for (const auto k : values) {
        const auto image = util::vec<T>(b[0] + i * box[0], b[1] + j * box[1],
                                        b[2] + k * box[2]);
        const auto diff =
            util::vec<T>(a[0] - image[0], a[1] - image[1], a[2] - image[2]);
        const auto dist2 = util::dot(diff, diff);
        if (dist2 < shortest) {
          closestImage = image;
          shortest = dist2;
        }
      }
    }
  }
  return closestImage;
}

std::vector<std::vector<util::rvec>> testCases() {
  auto cases = std::vector<std::vector<util::rvec>>(5);

  cases[0] = {
      {util::rvec(1, 0, 0)}, {util::rvec(8, 0, 0)}, {util::rvec(-2, 0, 0)}};
  cases[1] = {
      {util::rvec(1, 0, 0)}, {util::rvec(9, 0, 0)}, {util::rvec(-1, 0, 0)}};
  cases[2] = {
      {util::rvec(8, 0, 0)}, {util::rvec(1, 0, 0)}, {util::rvec(11, 0, 0)}};
  cases[3] = {
      {util::rvec(1, 1, 0)}, {util::rvec(8, 8, 0)}, {util::rvec(-2, -2, 0)}};
  cases[4] = {
      {util::rvec(1, 1, 1)}, {util::rvec(8, 8, 8)}, {util::rvec(-2, -2, -2)}};

  return cases;
}

const auto BOX = util::rvec(10, 10, 10);
const auto TESTCASES = testCases();

TEST(PBC, alternative_impl) {
  for (const auto test : TESTCASES) {
    const auto close = closestImageExhaustive(test[0], test[1], BOX);
    for (auto i = 0u; i < 3; ++i) {
      EXPECT_EQ(close[i], test[2][i]);
    }
  }
}

TEST(PBC, closest_image) {
  for (const auto test : TESTCASES) {
    const auto close = util::pbc::closestImage(test[0], test[1], BOX);
    for (auto i = 0u; i < 3; ++i) {
      EXPECT_EQ(close[i], test[2][i]);
    }
  }
}

TEST(PBC, closest_image_random_placements) {
  const auto n = 10000;
  auto rng = util::RNGEngine(RAND_SEED);
  const auto boxes = std::vector<util::rvec>{
      util::rvec(10, 10, 10), util::rvec(10, 20, 20), util::rvec(10, 20, 30)};

  for (const auto box : boxes) {
    for (auto i = 0; i < n; ++i) {
      const auto a = util::randomVec(box, rng);
      const auto b = util::randomVec(box, rng);
      const auto method_1 = closestImageExhaustive(a, b, box);
      const auto method_2 = util::pbc::closestImage(a, b, box);
      for (auto j = 0u; j < 3; ++j) {
        EXPECT_NEAR(method_1[j], method_2[j], 1e-14);
      }
    }
  }
}
