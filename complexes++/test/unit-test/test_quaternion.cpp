// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#include <cmath>
#include <random>

#include "complexes_test.h"
#include "util/array.h"
#include "util/moves.h"
#include "util/quaternions/quat.h"
#include "util/random.h"
#include "gtest/gtest.h"

namespace q = util::quaternions;

TEST(QUAT_TEST, norm) {
  auto Quater = 1. / std::sqrt(4);
  auto qa = q::Quaternion<double>(Quater, Quater, Quater, Quater);
  EXPECT_DOUBLE_EQ(1, qa.norm());
}

TEST(QUAT_TEST, conjugate) {
  auto qa = q::Quaternion<double>(0.5, 0.5, 0.5, 0.5).conjugate();
  EXPECT_EQ(q::Quaternion<double>(0.5, -0.5, -0.5, -0.5), qa);
}

TEST(QUAT_TEST, initialization_from_vector_and_angle) {
  auto const axis = util::rvec(1. / 3, 1. / 3, 1. / 3);
  auto q = q::Quaternion<double>(axis, 0);
  EXPECT_NEAR_VEC(q::rotateVector(q, axis), axis, 1e-8);

  q = q::Quaternion<double>(axis, 2 * M_PI);
  EXPECT_NEAR_VEC(q::rotateVector(q, axis), axis, 1e-8);

  q = q::Quaternion<double>(axis, M_PI);
  EXPECT_NEAR_VEC(q::rotateVector(q, axis), axis, 1e-8);

  auto const p = util::rvec(1.0, 0.0, 0.0);
  auto const axis2 = util::rvec(0.0, 0.0, 1.0);
  q = q::Quaternion<double>(axis2, M_PI / 2);
  EXPECT_NEAR_VEC(q::rotateVector(q, p), util::rvec(0.0, 1.0, 0.0), 1e-8);
}

TEST(QUAT_TEST, toMat) {
  auto const axis = util::rvec(1. / 3, 1. / 3, 1. / 3);
  auto q = q::Quaternion<double>(axis, 0).toMat();
  EXPECT_NEAR_VEC(util::rotateVector(q, axis), axis, 1e-8);

  q = q::Quaternion<double>(axis, 2 * M_PI).toMat();
  EXPECT_NEAR_VEC(util::rotateVector(q, axis), axis, 1e-8);

  q = q::Quaternion<double>(axis, M_PI).toMat();
  EXPECT_NEAR_VEC(util::rotateVector(q, axis), axis, 1e-8);

  auto const p = util::rvec(1.0, 0.0, 0.0);
  auto const axis2 = util::rvec(0.0, 0.0, 1.0);
  q = q::Quaternion<double>(axis2, M_PI / 2).toMat();
  EXPECT_NEAR_VEC(util::rotateVector(q, p), util::rvec(0.0, 1.0, 0.0), 1e-8);
}

TEST(QUAT_TEST, multiplication) {
  auto is2 = 1. / std::sqrt(2);
  auto right_angle = q::Quaternion<double>(is2, 0, 0, is2);
  auto half_turn = right_angle * right_angle;

  EXPECT_EQ(q::Quaternion<double>(0, 0, 0, 1), half_turn);

  auto full_turn = half_turn * half_turn;
  // Because of the double cover of SO(3) through the Quaternions we only need
  // to make sure that the magnitude is correct the sign doesn't matter.
  EXPECT_EQ(q::Quaternion<double>(0, 0, 0, 1), full_turn);
}

TEST(QUAT_TEST, unity) {
  const auto eye = q::unity<double>().toMat();
  for (auto i = 0; i < 3; ++i) {
    for (auto j = 0; j < 3; ++j) {
      if (i == j) {
        EXPECT_DOUBLE_EQ(1, eye(i, j));
      } else {
        EXPECT_DOUBLE_EQ(0, eye(i, j));
      }
    }
  }
}

util::rvec matrix_mul(const util::rArray& m, const util::rvec& v) {
  auto res = util::rvec(0, 0, 0);
  for (auto i = 0; i < 3; ++i) {
    for (auto j = 0; j < 3; ++j) {
      res[i] += m(i, j) * v[j];
    }
  }
  return res;
}

TEST(QUAT_TEST, rotation) {
  const auto is2 = 1. / std::sqrt(2);
  const auto right_angle = q::Quaternion<double>(is2, 0, 0, is2);
  const auto point = util::vec<double>(is2, is2, 0);

  const auto point2 = q::rotateVector(right_angle, point);

  EXPECT_DOUBLE_EQ(-is2, point2.at(0));
  EXPECT_DOUBLE_EQ(is2, point2.at(1));
  // this is a bit tricky. Floating-point comparison isn't easy and get really
  // complicated around 0. With normaliztion is gets really weird. To avoid all
  // of this we just add a constant offset.
  EXPECT_DOUBLE_EQ(1, point2.at(2) + 1);

  ////////////////////////////////////////////
  // Ok now to a more complicated rotation. //
  ////////////////////////////////////////////
  const auto angle = M_PI_2;
  const auto isr3 = 1. / std::sqrt(3);

  // quaternion
  const auto quat =
      q::Quaternion<double>(cos(angle * .5), sin(angle * .5) * isr3,
                            -sin(angle * .5) * isr3, sin(angle * .5) * isr3);

  // to further test it I will also calculate the rotation diffusion matrix for
  // this quaternion.
  auto mat = util::rArray(3, 3);
  const auto cosa = cos(angle);
  const auto sina = sin(angle);
  mat(0, 0) = cosa + (1 - cosa) / 3;
  mat(0, 1) = -(1 - cosa) / 3 - isr3 * sina;
  mat(0, 2) = (1 - cosa) / 3 - isr3 * sina;
  mat(1, 0) = -(1 - cosa) / 3 + isr3 * sina;
  mat(1, 1) = mat(0, 0);
  mat(1, 2) = -(1 - cosa) / 3 - isr3 * sina;
  mat(2, 0) = (1 - cosa) / 3 + isr3 * sina;
  mat(2, 1) = -(1 - cosa) / 3 + isr3 * sina;
  mat(2, 2) = mat(0, 0);

  // vector to rotate
  auto v = util::rvec(isr3, isr3, isr3);
  for (auto dummy = 0; dummy < 100; ++dummy) {
    const auto res = q::rotateVector(quat, v);
    const auto res_mat = matrix_mul(mat, v);
    v = res;
    EXPECT_NEAR_VEC(res_mat, res, .5e-10);
  }
}
