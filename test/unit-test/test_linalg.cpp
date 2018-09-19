#include "complexes_test.h"
#include "gtest/gtest.h"

#include "util/array.h"
#include "util/linalg.h"
#include "util/quaternions/quat.h"
#include "util/random.h"

TEST(LINALGTEST, matMul) {
  auto A = util::Matrix<double>();
  for (auto i = 0; i < 9; ++i) {
    A.data()[i] = i + 1;
  }
  auto res = util::Matrix<double>();
  const auto res_vals =
      std::array<double, 9>{{30, 36, 42, 66, 81, 96, 102, 126, 150}};
  for (auto i = 0; i < 9; ++i) {
    res.data()[i] = res_vals[i];
  }

  const auto C = util::matMul(A, A);

  for (auto i = 0; i < 3; ++i) {
    for (auto j = 0; j < 3; ++j) {
      EXPECT_DOUBLE_EQ(res(i, j), C(i, j))
          << fmt::format(" at i={}, j={}", i, j);
    }
  }
}

template <typename RealType>
util::Matrix<RealType> random3x3(util::RNGEngine& rng) {
  auto dist = std::uniform_real_distribution<RealType>{0, 1};
  auto arr = util::Matrix<RealType>();
  for (auto i = 0; i < 3; ++i) {
    for (auto j = 0; j < 3; ++j) {
      arr(i, j) = dist(rng);
    }
  }
  return arr;
}

TEST(LINALGTEST, inv) {
  auto rng = util::RNGEngine(42424242424242);
  const auto maxTrials = 1000;

  // See if inverse times transpose equal the identity matrix
  for (auto k = 0; k < maxTrials; ++k) {
    const auto A = random3x3<double>(rng);
    const auto Ai = util::inv(A);
    const auto eye = util::matMul(A, Ai);

    for (auto i = 0; i < 3; ++i) {
      for (auto j = 0; j < 3; ++j) {
        if (i == j) {
          EXPECT_NEAR(eye(i, j), 1, 1e-10)
              << fmt::format(" at i={}, j={}", i, j);
        } else {
          EXPECT_NEAR(eye(i, j), 0, 1e-10)
              << fmt::format(" at i={}, j={}", i, j);
        }
      }
    }
  }

  // for rotation matrices the inverse equals the transpose
  for (auto k = 0; k < maxTrials; ++k) {
    const auto A = util::randomRotation(3.0, rng);
    const auto Ai = util::inv(A);
    const auto At = util::transpose(A);

    for (auto i = 0; i < 3; ++i) {
      for (auto j = 0; j < 3; ++j) {
        EXPECT_NEAR(Ai(i, j), At(i, j), 1e-10)
            << fmt::format(" at i={}, j={}", i, j);
      }
    }
  }
}

TEST(LINALGTEST, det) {
  auto rng = util::RNGEngine(42424242424242);
  const auto maxTrials = 1000;
  for (auto k = 0; k < maxTrials; ++k) {
    EXPECT_NEAR(util::det(util::randomRotation(3.0, rng)), 1, 1e-10);
  }
}

TEST(LINALGTEST, transpose) {
  auto rng = util::RNGEngine(42424242424242);
  const auto A = random3x3<double>(rng);
  const auto At = util::transpose(A);
  for (auto i = 0; i < 3; ++i) {
    for (auto j = 0; j < 3; ++j) {
      EXPECT_DOUBLE_EQ(A(i, j), At(j, i));
    }
  }
}

TEST(LINALGTEST, scale) {
  const auto box = util::rvec(1, 2, 3);
  const auto s = 2.0;
  const auto newbox = util::scale(box, s);
  EXPECT_EQ_VEC(newbox, util::rvec(2, 4, 6));

  auto ar = util::Array<double>(2, 3);
  for (auto i = 0; i < ar.rows(); ++i) {
    for (auto j = 0; j < ar.cols(); ++j) {
      ar(i, j) = i + j;
    }
  }

  const auto new_ar = util::scale(ar, s);
  for (auto i = 0; i < ar.rows(); ++i) {
    for (auto j = 0; j < ar.cols(); ++j) {
      EXPECT_DOUBLE_EQ(new_ar(i, j), s * (i + j));
    }
  }
}

TEST(LINALGTEST, volume) {
  const auto box = util::rvec(2, 3, 4);
  EXPECT_DOUBLE_EQ(util::volume(box), 2 * 3 * 4);
}
