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
#ifndef LINALG_H
#define LINALG_H

#include <numeric>

#include "util/array.h"
#include "util/util.h"

namespace util {

// attention, this function does not really return the distance
// but the vector that is the difference between the two input vectors.
void distance(util::rvec &d, const int i, const int j, const util::rArray &xyz);
// returns the distance between vectors v1 and v2
double distance(const util::rvec &v1, const util::rvec &v2);
void normalize(util::rvec &v);
rvec normalizeVector(const rvec &xyz);

template <typename Real,
          typename = std::enable_if<std::is_floating_point<Real>::value>>
Matrix<Real> matMul(const Matrix<Real> &a, const Matrix<Real> &b) {
  auto res = Matrix<Real>(0);
  for (auto i = 0; i < 3; ++i) {
    for (auto j = 0; j < 3; ++j) {
      for (auto k = 0; k < 3; ++k) {
        res(i, j) += a(i, k) * b(k, j);
      }
    }
  }
  return res;
}

template <typename Real, int rows, int cols,
          typename = std::enable_if<std::is_floating_point<Real>::value>>
FixedArray<Real, cols, rows> transpose(const FixedArray<Real, rows, cols> &a) {
  auto b = FixedArray<Real, cols, rows>();
  for (auto i = 0; i < cols; ++i) {
    for (auto j = 0; j < rows; ++j) {
      b(i, j) = a(j, i);
    }
  }
  return b;
}

template <typename Real,
          typename = std::enable_if<std::is_floating_point<Real>::value>>
Real dot(const vec<Real> &v1, const vec<Real> &v2) {
  return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

template <typename Real,
          typename = std::enable_if<std::is_floating_point<Real>::value>>
void cross(vec<Real> &out, const vec<Real> &v1, const vec<Real> &v2) {
  out[0] = v1[1] * v2[2] - v1[2] * v2[1];
  out[1] = v1[2] * v2[0] - v1[0] * v2[2];
  out[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

template <typename Real,
          typename = std::enable_if<std::is_floating_point<Real>::value>>
Real _det(const vec<Real> &x, const vec<Real> &y, const vec<Real> &z) {
  auto c = vec<Real>(0, 0, 0);
  cross(c, y, z);
  return dot(x, c);
}

// determinante for 3x3 matrix
template <typename Real,
          typename = std::enable_if<std::is_floating_point<Real>::value>>
Real det(const util::Matrix<Real> &A) {
  const auto x = vec<Real>(A(0, 0), A(1, 0), A(2, 0));
  const auto y = vec<Real>(A(0, 1), A(1, 1), A(2, 1));
  const auto z = vec<Real>(A(0, 2), A(1, 2), A(2, 2));
  return _det(x, y, z);
}

// inverse of 3x3 matrix
template <typename Real,
          typename = std::enable_if<std::is_floating_point<Real>::value>>
util::Matrix<Real> inv(const util::Matrix<Real> &A) {
  const auto x = vec<Real>(A(0, 0), A(1, 0), A(2, 0));
  const auto y = vec<Real>(A(0, 1), A(1, 1), A(2, 1));
  const auto z = vec<Real>(A(0, 2), A(1, 2), A(2, 2));

  auto res = Matrix<Real>();

  auto c = vec<Real>(0, 0, 0);
  cross(c, y, z);
  res(0, 0) = c[0];
  res(0, 1) = c[1];
  res(0, 2) = c[2];
  cross(c, z, x);
  res(1, 0) = c[0];
  res(1, 1) = c[1];
  res(1, 2) = c[2];
  cross(c, x, y);
  res(2, 0) = c[0];
  res(2, 1) = c[1];
  res(2, 2) = c[2];

  const auto d = _det(x, y, z);

  for (auto i = 0; i < 3; ++i) {
    for (auto j = 0; j < 3; ++j) {
      res(i, j) /= d;
    }
  }

  return res;
}

template <typename Number,
          typename = std::enable_if<std::is_arithmetic<Number>::value>>
std::vector<Number> arange(const int n) {
  std::vector<Number> arr(n);
  std::iota(std::begin(arr), std::end(arr), 0);
  return arr;
}

// I want to sort one array and retain the indices in order to be able
// to draw things from a second array using that ordering.
// example usage:
// for (auto i: sortIndices(v)) {
//  cout << v[i] << endl;
//}
template <typename T> std::vector<size_t> sortIndices(const std::vector<T> &v) {
  auto idx = arange<std::size_t>(v.size());
  // sort indexes based on comparing values in v
  sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) { return v[i1] < v[i2]; });

  return idx;
}

// TODO: can I compile time introspect the type T contained in the vector?
template <typename T> std::vector<T> diff(const std::vector<T> &v) {
  DEBUG_ASSERT(v.size() >= 2,
               "trying to diff a vector with one element or less.");
  auto d = std::vector<T>(v.size() - 1);
  for (auto i = 0u; i < v.size() - 1; ++i) {
    d[i] = v[i + 1] - v[i];
  }
  return d;
}

// TODO: allow arbitrary container with elements of type int
template <typename T> bool isConsecutive(const std::vector<T> &c) {
  DEBUG_ASSERT(c.size() > 0, "Can't determine if 0 elements are consecutive");
  if (c.size() == 1) {
    return true;
  }

  const auto difference = diff(c);
  for (const auto el : difference) {
    if (el != 1) {
      return false;
    }
  }

  return true;
}

template <class ContainerType>
auto upperTriangleSum(const ContainerType &a) noexcept ->
    typename std::enable_if<
        ContainerType::IsRowMajor,
        typename std::remove_reference<decltype(a(0, 0))>::type>::type {
  auto sum = decltype(a(0, 0)){0};
  for (auto j = 0; j < a.cols(); ++j) {
    for (auto i = j; i < a.rows(); ++i) {
      sum += a(i, j);
    }
  }
  return sum;
}

template <class ContainerType>
auto upperTriangleSum(const ContainerType &a) noexcept ->
    typename std::enable_if<
        !ContainerType::IsRowMajor,
        typename std::remove_reference<decltype(a(0, 0))>::type>::type {
  auto sum = decltype(a(0, 0)){0};
  for (auto i = 0; i < a.rows(); ++i) {
    for (auto j = i; j < a.cols(); ++j) {
      sum += a(i, j);
    }
  }
  return sum;
}

template <typename Real,
          typename = std::enable_if<std::is_floating_point<Real>::value>>
Real volume(const util::vec<Real> &box) {
  return box[0] * box[1] * box[2];
}

template <typename Real,
          typename = std::enable_if<std::is_floating_point<Real>::value>>
vec<Real> scale(const vec<Real> &v, const Real fac) {
  return util::vec<Real>(v[0] * fac, v[1] * fac, v[2] * fac);
}

template <typename Real,
          typename = std::enable_if<std::is_floating_point<Real>::value>>
Array<Real> scale(const Array<Real> &a, const Real fac) {
  const auto rows = a.rows();
  const auto cols = a.cols();
  auto b = util::Array<Real>(rows, cols);
  for (auto i = 0; i < rows; ++i) {
    for (auto j = 0; j < cols; ++j) {
      b(i, j) = fac * a(i, j);
    }
  }
  return b;
}

} // namespace util

#endif // LINALG_H
