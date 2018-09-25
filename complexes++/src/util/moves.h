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
#ifndef MOVES_H
#define MOVES_H

#include "util/array.h"
#include "util/linalg.h"
#include "util/quaternions/quat.h"

namespace util {

template <typename Real,
          typename = std::enable_if<std::is_floating_point<Real>::value>>
void translateByConstantInPlace(const Real constant, Array<Real> *xyzToShift) {
  for (auto i = 0; i < (*xyzToShift).rows(); ++i) {
    for (auto j = 0; j < 3; ++j) {
      (*xyzToShift)(i, j) += constant;
    }
  }
}

template <typename Real,
          typename = std::enable_if<std::is_floating_point<Real>::value>>
void translateByConstantInPlace(const vec<Real> &constants,
                                Array<Real> *xyzToShift) {
  for (auto i = 0; i < (*xyzToShift).rows(); ++i) {
    for (auto j = 0; j < 3; ++j) {
      (*xyzToShift)(i, j) += constants[j];
    }
  }
}

template <typename Real,
          typename = std::enable_if<std::is_floating_point<Real>::value>>
Array<Real> translateByConstant(const Array<Real> &xyz,
                                const vec<Real> &constant) {
  auto res = xyz;
  translateByConstantInPlace(constant, &res);
  return res;
}

template <typename Real,
          typename = std::enable_if<std::is_floating_point<Real>::value>>
Array<Real> rotateWithQuat(const Array<Real> &xyz, const vec<Real> &com,
                           const quaternions::Quaternion<Real> &quat) {
  auto res = Array<Real>(xyz.rows(), xyz.cols());
  for (auto i = 0; i < xyz.rows(); ++i) {
    auto vec = rvec(xyz(i, 0) - com[0], xyz(i, 1) - com[1], xyz(i, 2) - com[2]);
    vec = quaternions::rotateVector(quat, vec);
    for (auto j = 0; j < 3; ++j) {
      res(i, j) = vec[j] + com[j];
    }
  }
  return res;
}

template <typename Real,
          typename = std::enable_if<std::is_floating_point<Real>::value>>
vec<Real> rotateVector(const Matrix<Real> &R, const vec<Real> &v) {
  auto u = vec<Real>(0, 0, 0);
  for (auto j = 0; j < 3; ++j) {
    for (auto k = 0; k < 3; ++k) {
      u[j] += R(j, k) * v[k];
    }
  }
  return u;
}

template <typename Real,
          typename = std::enable_if<std::is_floating_point<Real>::value>>
void rotateWithMatInPlace(const vec<Real> &com, const Matrix<Real> &mat,
                          Array<Real> *xyz) {
  for (auto i = 0; i < (*xyz).rows(); ++i) {
    const auto v = rotateVector(mat, vec<Real>((*xyz)(i, 0) - com[0],
                                               (*xyz)(i, 1) - com[1],
                                               (*xyz)(i, 2) - com[2]));
    for (auto j = 0; j < 3; ++j) {
      (*xyz)(i, j) = v[j] + com[j];
    }
  }
}

template <typename Real,
          typename = std::enable_if<std::is_floating_point<Real>::value>>
Array<Real> rotateWithMat(const Array<Real> &xyz, const vec<Real> &com,
                          const Matrix<Real> &mat) {
  auto res = xyz;
  rotateWithMatInPlace(com, mat, &res);
  return res;
}

template <typename Real,
          typename = std::enable_if<std::is_floating_point<Real>::value>>
vec<Real> centroid(const Array<Real> &xyz) {
  auto com = vec<Real>(0);
  for (auto i = 0; i < xyz.rows(); ++i) {
    for (auto j = 0; j < 3; ++j) {
      com[j] += xyz(i, j);
    }
  }
  for (auto j = 0; j < 3; ++j) {
    com[j] /= xyz.rows();
  }
  return com;
}

template <typename Real,
          typename = std::enable_if<std::is_floating_point<Real>::value>>
vec<Real> distanceVec(const vec<Real> &v1, const vec<Real> &v2) {
  return util::vec<Real>(v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2]);
}

template <typename Real,
          typename = std::enable_if<std::is_floating_point<Real>::value>>
Array<Real> scaleCentroid(const Array<Real> &xyz, const Real fac) {
  const auto cen = centroid(xyz);
  const auto newCen = scale(cen, fac);
  const auto transVec = distanceVec(newCen, cen);
  return translateByConstant(xyz, transVec);
}
} // namespace util

#endif // MOVES_H
