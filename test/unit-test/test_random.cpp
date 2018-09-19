// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#include <cmath>

#include "complexes_test.h"
#include "gtest/gtest.h"

#include "util/linalg.h"
#include "util/random.h"

TEST(RANDOM, state) {
  auto rng = util::RNGEngine(RAND_SEED);
  for (auto i = 0; i < 100; ++i) {
    rng();
  }
  const auto old_state = util::getRNGState(rng);
  const auto next = rng();
  util::setRNGState(rng, old_state);
  const auto other_next = rng();

  EXPECT_EQ(next, other_next);
}

TEST(RANDOM, indices_random_order) {
  auto rng = util::RNGEngine(RAND_SEED);
  const auto N = 1000u;
  const auto rand = util::randomOrderIndices(N, rng);

  auto order = true;
  for (auto i = 0u; i < N; ++i) {
    order = order && (i == rand[i]);
  }

  EXPECT_FALSE(order);
}

TEST(RANDOM, rand_vec) {
  auto rng = util::RNGEngine(RAND_SEED);
  const auto scale = util::rvec(10, 20, 30);
  for (auto j = 0; j < 100; ++j) {
    const auto vec = util::randomVec(scale, rng);
    for (auto i = 0u; i < 3; ++i) {
      EXPECT_TRUE(vec[i] < scale[i] && vec[i] > -scale[i]);
    }
  }
}

TEST(RANDOM, random_quat) {
  const auto v = util::vec<double>(.1, .1, .1);
  auto rng = util::RNGEngine(util::fullEntropySeed(42));

  for (auto i = 0; i < 10000; ++i) {
    const auto qa =
        util::randomQuat(v, util::rvec(.0011, .0011, .0011), 1.0, rng);
    EXPECT_DOUBLE_EQ(1, qa.norm());
  }
}

template <typename T>
class UniformCDF {
 public:
  UniformCDF(const T a, const T b) : m_a(a), m_b(b) {}
  T operator()(const T x) const { return (x - m_a) / (m_b - m_a); }

 private:
  const T m_a;
  const T m_b;
};

TEST(RANDOM, rand_vec_kstest) {
  auto rng = util::RNGEngine(RAND_SEED);
  const auto scale = util::rvec(10, 20, 30);
  const auto N = 10000;
  auto rvs = std::array<std::vector<double>, 3>();
  for (auto& vec : rvs) {
    vec.reserve(N);
  }

  for (auto j = 0; j < N; ++j) {
    const auto vec = util::randomVec(scale, rng);
    for (auto i = 0u; i < rvs.size(); ++i) {
      rvs[i].push_back(vec[i]);
    }
  }

  for (auto i = 0u; i < rvs.size(); ++i) {
    const auto p_value =
        kstest<double>(rvs[i], UniformCDF<double>(-scale[i], scale[i]));
    // 50 % of the time we expect the difference between the CDF estimated from
    // rvs and the analytic cdf is as large as the observed one this time. This
    // mean the NULL-hypothesis is confirmed and we really sample a normal
    // distribution.
    EXPECT_TRUE(p_value > .4);
    // For detail see
    // https://sites.google.com/a/ucsc.edu/krumholz/teaching-and-courses/ast119_w15/class-10
  }
}
