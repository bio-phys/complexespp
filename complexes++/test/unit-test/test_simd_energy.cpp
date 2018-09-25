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
#include "complexes_test.h"
#include "gtest/gtest.h"
#include <cmath>

#include "constants.h"
#include "util/array.h"
#include "vectorization/functionsvec.h"

#include "InastempConfig.h"
#include "SCALAR/InaVecSCALARDouble.hpp"

#ifdef INASTEMP_USE_AVX2
#include "AVX2/InaVecAVX2Double.hpp"
#endif

#ifdef INASTEMP_USE_SSE3
#include "SSE3/InaVecSSE3Double.hpp"
#endif

namespace nc = constants::natural;

// NOTE ABOUT THE UNITS:
// length [Anstroem]
// energy [kT]

#define SIMD_EXPECT_DOUBLE_EQ(VecTypeClass, X, Y)                              \
  {                                                                            \
    const VecTypeClass v1 = X;                                                 \
    const VecTypeClass v2 = Y;                                                 \
    for (int idx = 0; idx < VecType::VecLength; ++idx) {                       \
      EXPECT_DOUBLE_EQ(v1.at(idx), v2.at(idx));                                \
    }                                                                          \
  }

#define SIMD_EXPECT_NEAR(VecTypeClass, X, Y, Z)                                \
  {                                                                            \
    const VecTypeClass v1 = X;                                                 \
    const VecTypeClass v2 = Y;                                                 \
    for (int idx = 0; idx < VecTypeClass::VecLength; ++idx) {                  \
      EXPECT_NEAR(v1.at(idx), v2.at(idx), Z);                                  \
    }                                                                          \
  }

#define SIMD_EXPECT_TRUE(X) EXPECT_TRUE((X).isAllTrue());

// Macro to test SIMD test function for several types
#if defined INASTEMP_USE_AVX2
#define TEST_SIMD_FUNC(name, func)                                             \
  TEST(SIMDENERGYTEST, name) {                                                 \
    func<InaVecSCALAR<double>>();                                              \
    func<InaVecBestType<double>>();                                            \
    func<InaVecAVX2<double>>();                                                \
  }
#elif defined INASTEMP_USE_SSE3
#define TEST_SIMD_FUNC(name, func)                                             \
  TEST(SIMDENERGYTEST, name) {                                                 \
    func<InaVecSCALAR<double>>();                                              \
    func<InaVecBestType<double>>();                                            \
    func<InaVecSSE3<double>>();                                                \
  }
#else
#define TEST_SIMD_FUNC(name, func)                                             \
  TEST(SIMDENERGYTEST, name) {                                                 \
    func<InaVecSCALAR<double>>();                                              \
    func<InaVecBestType<double>>();                                            \
  }
#endif

template <class VecType> void simd_lennard_jones() {
  const VecType r0 = std::pow(2, 1. / 6.);
  const VecType sigma = 1.;

  // There isn't much to test the implementation with out writing the function
  // twice. So lets just check that the extrema and other characterics are OK.

  VecType epsilon = -2;
  SIMD_EXPECT_DOUBLE_EQ(
      VecType, epsilon,
      simd::util::lennardJones<VecType>(r0 * r0, epsilon, sigma));
  SIMD_EXPECT_TRUE(simd::util::lennardJones<VecType>((2. * r0) * (2. * r0),
                                                     epsilon, sigma) < 0);
  SIMD_EXPECT_TRUE(simd::util::lennardJones<VecType>(
                       VecType(.99) * VecType(.99), epsilon, sigma) > 0);
  SIMD_EXPECT_DOUBLE_EQ(
      VecType, 0.,
      simd::util::lennardJones<VecType>(VecType(1.), epsilon, sigma));

  epsilon = 2.;
  SIMD_EXPECT_DOUBLE_EQ(
      VecType, epsilon,
      simd::util::lennardJones<VecType>(r0 * r0, epsilon, sigma));
  SIMD_EXPECT_TRUE(simd::util::lennardJones<VecType>(r0 * r0, epsilon, sigma) >
                   0);
  SIMD_EXPECT_TRUE(simd::util::lennardJones<VecType>(VecType(std::pow(3, 2)),
                                                     epsilon, sigma) > 0);
  SIMD_EXPECT_DOUBLE_EQ(
      VecType, 2 * epsilon,
      simd::util::lennardJones<VecType>(VecType(1.), epsilon, sigma));

  epsilon = 0.;
  SIMD_EXPECT_DOUBLE_EQ(
      VecType, VecType(.01 / std::pow(std::pow(2, 1. / 6.), 12)),
      simd::util::lennardJones<VecType>(r0 * r0, epsilon, sigma));
  SIMD_EXPECT_TRUE(simd::util::lennardJones<VecType>(r0 * r0, epsilon, sigma) >
                   0);
}
TEST_SIMD_FUNC(lennard_jones, simd_lennard_jones)

template <class VecType> void simd_debye_hueckel() {
  const VecType charge1 = 1;
  const VecType r = 1;
  const VecType debyeLength = 1;
  const VecType D = nc::elementaryCharge * nc::elementaryCharge /
                    (4 * M_PI * nc::epsilon_0 * constants::units::energy *
                     constants::units::angstrom);

  SIMD_EXPECT_TRUE(
      simd::util::debyeHueckel<VecType>(charge1, 1, r, debyeLength, D) > 0);
  SIMD_EXPECT_TRUE(
      simd::util::debyeHueckel<VecType>(charge1, -1, r, debyeLength, D) < 0);
  SIMD_EXPECT_NEAR(
      VecType, 0.,
      simd::util::debyeHueckel<VecType>(charge1, 0., r, debyeLength, D), 1E-12);
  SIMD_EXPECT_NEAR(
      VecType, exp(-1),
      simd::util::debyeHueckel<VecType>(charge1, 1, r, debyeLength, D), 1E-12);
}
TEST_SIMD_FUNC(debye_hueckel, simd_debye_hueckel)

template <class VecType> void simd_smooth() {
  const VecType rmin_pow2 = 1.5 * 1.5;
  const VecType rmax_pow2 = 2 * 2;
  const VecType sigma = 1;
  SIMD_EXPECT_NEAR(
      VecType, 0,
      simd::util::smooth<VecType>(VecType(5), sigma, rmin_pow2, rmax_pow2),
      1E-12);
  SIMD_EXPECT_NEAR(
      VecType, 1,
      simd::util::smooth<VecType>(VecType(1), sigma, rmin_pow2, rmax_pow2),
      1E-12);
  SIMD_EXPECT_NEAR(VecType, .553480,
                   simd::util::smooth<VecType>(VecType(1.75 * 1.75), sigma,
                                               rmin_pow2, rmax_pow2),
                   1E-5);
}
TEST_SIMD_FUNC(smooth, simd_smooth)
