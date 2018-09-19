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
#include "util/moves.h"
#include "util/pbc.h"

namespace q = util::quaternions;

TEST(MOVESTEST, applyBPC) {
  auto xyz = util::rArray(1, 3);
  const auto box = util::rvec(5, 5, 5);

  xyz(0, 0) = 4;
  xyz(0, 1) = 6;
  xyz(0, 2) = -1;

  const auto xyz_pbc = util::pbc::applyPBC(xyz, box);

  EXPECT_DOUBLE_EQ(4, xyz_pbc(0, 0));
  EXPECT_DOUBLE_EQ(1, xyz_pbc(0, 1));
  EXPECT_DOUBLE_EQ(4, xyz_pbc(0, 2));
}

TEST(MOVESTEST, ToClosest) {
  auto d = util::rvec(6, 1, -8);
  const auto box = util::rvec(10, 10, 10);

  util::pbc::ToClosest(&d, box);

  EXPECT_DOUBLE_EQ(-4, d.at(0));
  EXPECT_DOUBLE_EQ(1, d.at(1));
  EXPECT_DOUBLE_EQ(2, d.at(2));

  const auto box2 = util::rvec(20, 20, 20);
  d = util::rvec(5, 5, 5);
  util::pbc::ToClosest(&d, box2);

  EXPECT_DOUBLE_EQ(5, d.at(0));
  EXPECT_DOUBLE_EQ(5, d.at(1));
  EXPECT_DOUBLE_EQ(5, d.at(2));
}

TEST(MOVES_TEST, translateByConstant) {
  auto xyz = util::rArray(1, 3);
  const auto constant = util::rvec(5, 5, 5);

  xyz(0, 0) = 4;
  xyz(0, 1) = 6;
  xyz(0, 2) = -1;

  const auto trans = util::translateByConstant(xyz, constant);

  EXPECT_DOUBLE_EQ(9, trans(0, 0));
  EXPECT_DOUBLE_EQ(11, trans(0, 1));
  EXPECT_DOUBLE_EQ(4, trans(0, 2));
}

TEST(MOVES_TEST, rotateWithQuat) {
  auto angle = M_PI_2;
  auto isr3 = 1. / std::sqrt(3);

  // quaternion
  auto quat =
      q::Quaternion<double>(cos(angle * .5), sin(angle * .5) * isr3,
                            -sin(angle * .5) * isr3, sin(angle * .5) * isr3);

  const auto v_org = util::rvec(isr3, isr3, isr3);
  const auto v_org2 = util::rvec(isr3, isr3, -isr3);

  // Generate system v_org and v_org2 as coordinates
  auto xyz = util::rArray(2, 3);

  xyz(0, 0) = v_org[0];
  xyz(0, 1) = v_org[1];
  xyz(0, 2) = v_org[2];
  xyz(1, 0) = v_org2[0];
  xyz(1, 1) = v_org2[1];
  xyz(1, 2) = v_org2[2];

  // shift by a constant com;
  const auto com = util::rvec(2, 2, 2);
  for (auto i = 0; i < xyz.rows(); ++i) {
    for (auto j = 0; j < xyz.cols(); ++j) {
      xyz(i, j) += com.at(j);
    }
  }

  const auto erg = util::rotateWithQuat(xyz, com, quat);
  const auto v_res1 = q::rotateVector(quat, v_org);
  const auto v_res2 = q::rotateVector(quat, v_org2);

  EXPECT_DOUBLE_EQ(v_res1[0] + com[0], erg(0, 0));
  EXPECT_DOUBLE_EQ(v_res1[1] + com[1], erg(0, 1));
  EXPECT_DOUBLE_EQ(v_res1[2] + com[2], erg(0, 2));

  EXPECT_DOUBLE_EQ(v_res2[0] + com[0], erg(1, 0));
  EXPECT_DOUBLE_EQ(v_res2[1] + com[1], erg(1, 1));
  EXPECT_DOUBLE_EQ(v_res2[2] + com[2], erg(1, 2));
}

TEST(MOVES_TEST, rotateWithMat) {
  auto angle = M_PI_2;
  auto isr3 = 1. / sqrt(3);

  // quaternion
  auto mat =
      q::Quaternion<double>(cos(angle * .5), sin(angle * .5) * isr3,
                            -sin(angle * .5) * isr3, sin(angle * .5) * isr3)
          .toMat();

  const auto v_org = util::rvec(isr3, isr3, isr3);
  const auto v_org2 = util::rvec(isr3, isr3, -isr3);

  // Generate system v_org and v_org2 as coordinates
  auto xyz = util::rArray(2, 3);

  xyz(0, 0) = v_org[0];
  xyz(0, 1) = v_org[1];
  xyz(0, 2) = v_org[2];
  xyz(1, 0) = v_org2[0];
  xyz(1, 1) = v_org2[1];
  xyz(1, 2) = v_org2[2];

  // shift by a constant com;
  const auto com = util::rvec(2, 2, 2);
  for (auto i = 0; i < xyz.rows(); ++i) {
    for (auto j = 0; j < xyz.cols(); ++j) {
      xyz(i, j) += com.at(j);
    }
  }

  const auto erg = util::rotateWithMat(xyz, com, mat);
  const auto v_res1 = util::rotateVector(mat, v_org);
  const auto v_res2 = util::rotateVector(mat, v_org2);

  EXPECT_DOUBLE_EQ(v_res1[0] + com[0], erg(0, 0));
  EXPECT_DOUBLE_EQ(v_res1[1] + com[1], erg(0, 1));
  EXPECT_DOUBLE_EQ(v_res1[2] + com[2], erg(0, 2));

  EXPECT_DOUBLE_EQ(v_res2[0] + com[0], erg(1, 0));
  EXPECT_DOUBLE_EQ(v_res2[1] + com[1], erg(1, 1));
  EXPECT_DOUBLE_EQ(v_res2[2] + com[2], erg(1, 2));
}

TEST(MOVES_TEST, rotateVector) {
  const auto eye = q::unity<double>().toMat();
  const auto v1 = util::vec<double>(1, 1, 1);
  const auto v2 = util::rotateVector(eye, v1);
  EXPECT_EQ_VEC(v1, v2);

  const auto is2 = 1. / std::sqrt(2);
  const auto right_angle = q::Quaternion<double>(is2, 0, 0, is2).toMat();
  const auto point = util::vec<double>(is2, is2, 0);
  const auto point2 = util::rotateVector(right_angle, point);

  EXPECT_DOUBLE_EQ(-is2, point2.at(0));
  EXPECT_DOUBLE_EQ(is2, point2.at(1));
  // this is a bit tricky. Floating-point comparison isn't easy and get really
  // complicated around 0. With normaliztion it gets really weird. To avoid all
  // of this we just add a constant offset.
  EXPECT_DOUBLE_EQ(1, point2.at(2) + 1);

  ////////////////////////////////////////////
  // Ok now to a more complicated rotation. //
  ////////////////////////////////////////////
  const auto angle = M_PI_2;
  const auto isr3 = 1. / std::sqrt(3);
  const auto quat =
      q::Quaternion<double>(cos(angle * .5), sin(angle * .5) * isr3,
                            -sin(angle * .5) * isr3, sin(angle * .5) * isr3);
  const auto mat = quat.toMat();

  // vector to rotate
  auto v = util::rvec(isr3, isr3, isr3);
  for (auto k = 0; k < 100; ++k) {
    const auto res = q::rotateVector(quat, v);
    const auto res_mat = util::rotateVector(mat, v);
    v = res;
    EXPECT_NEAR_VEC(res_mat, res, .5e-10);
  }

  // use known rotation to check against
  auto M = util::Matrix<double>();
  M(0, 0) = 1;
  M(1, 2) = -1;
  M(2, 1) = 1;
  v = util::rvec(1, 2, 3);
  const auto res = util::rotateVector(M, v);
  EXPECT_NEAR_VEC(res, util::rvec(1, -3, 2), .5e-10);
}

TEST(MOVES_TEST, centroid) {
  auto xyz = util::rArray(2, 3);

  xyz(0, 0) = 10;
  xyz(1, 0) = 20;

  xyz(0, 1) = 10;
  xyz(1, 1) = -10;

  xyz(0, 2) = 20;
  xyz(1, 2) = -10;

  const auto com = util::centroid(xyz);
  EXPECT_DOUBLE_EQ(15, com[0]);
  EXPECT_DOUBLE_EQ(0, com[1]);
  EXPECT_DOUBLE_EQ(5, com[2]);
}

TEST(TEST_MOVES, distance) {
  const auto a = util::rvec(1, 2, 3);
  const auto b = util::rvec(2, 3, 4);

  EXPECT_EQ_VEC(util::distanceVec(a, b), util::rvec(-1, -1, -1));
  EXPECT_EQ_VEC(util::distanceVec(b, a), util::rvec(1, 1, 1));
  EXPECT_EQ_VEC(util::distanceVec(b, b), util::rvec(0, 0, 0));
}

TEST(TEST_MOVES, scaleCentroid) {
  auto ar = util::Array<double>(2, 3);
  for (auto i = 0; i < ar.rows(); ++i) {
    for (auto j = 0; j < ar.cols(); ++j) {
      ar(i, j) = i * 2;
    }
  }
  auto s = 2.0;
  auto newar = util::scaleCentroid(std::move(ar), s);
  const auto center = util::centroid(newar);

  EXPECT_NEAR_VEC(center, util::rvec(1, 1, 1), 12);
}
