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
#ifndef ARRAY_H
#define ARRAY_H

#include <fmt/format.h>
#include <fmt/ostream.h>
#include <initializer_list>
#include <iostream>
#include <vector>

#include "io/serializer.h"
#include "util/util.h"

namespace util {

template <class T> class vec : public io::AbstractSerializable {
private:
  std::array<T, 3> m_data;

public:
  vec() : m_data{T(), T(), T()} {}
  vec(const T t) : m_data{t, t, t} {}
  vec(const T a, const T b, const T c) : m_data{a, b, c} {}

  // support moving
  vec(vec &&rhs) = default;
  vec &operator=(vec &&rhs) = default;
  // support copying
  vec(const vec &rhs) = default;
  vec &operator=(const vec &rhs) = default;

  T &operator[](int i) { return m_data[i]; }
  const T &operator[](int i) const { return m_data[i]; }
  const T &at(int i) const { return m_data.at(i); }

  int size() const { return 3; }

  vec &operator+=(const vec &inVec) {
    for (int idx = 0; idx < 3; ++idx) {
      m_data[idx] += inVec[idx];
    }
    return (*this);
  }

  vec &operator-=(const vec &inVec) {
    for (int idx = 0; idx < 3; ++idx) {
      m_data[idx] -= inVec[idx];
    }
    return (*this);
  }

  /// Serialization
  void serialize(io::Serializer &serializer) const final {
    serializer.append(m_data, "m_data");
  }
  vec(io::Deserializer &deserializer)
      : m_data(deserializer.restore<decltype(m_data)>("m_data")) {}
};
using rvec = vec<double>;

class RowMajor {
public:
  static const bool IsRowMajor = true;

  inline static int IndexOf(const int idxRow, const int idxCol,
                            const int inNbRows, const int inNbCols) {
    UNUSED(inNbCols);
    return inNbRows * idxCol + idxRow;
  }
};

class ColumnMajor {
public:
  static const bool IsRowMajor = false;

  inline static int IndexOf(const int idxRow, const int idxCol,
                            const int inNbRows, const int inNbCols) {
    UNUSED(inNbRows);
    return inNbCols * idxRow + idxCol;
  }
};

/* The Array class stores the values contiguously (all rows/cols in one block).
 */
template <typename T, typename Ordering = RowMajor>
class Array : public io::AbstractSerializable {
public:
  static const bool IsRowMajor = Ordering::IsRowMajor;

  explicit Array(int rows_, int cols_, const T &defaultValue = T())
      : m_rows(rows_), m_cols(cols_), m_pData(m_rows * m_cols, defaultValue) {}
  explicit Array(std::size_t rows_, std::size_t cols_,
                 const T &defaultValue = T())
      : m_rows(static_cast<int>(rows_)), m_cols(static_cast<int>(cols_)),
        m_pData(m_rows * m_cols, defaultValue) {}
  // support moving
  Array(Array &&rhs) = default;
  Array &operator=(Array &&rhs) = default;
  // support copying
  Array(const Array &rhs) = default;
  Array &operator=(const Array &rhs) = default;

  T &operator()(int row, int col) {
    DEBUG_ASSERT(row < m_rows && col < m_cols,
                 "Out of bound array access ({} {}) with cols={}, rows={}", row,
                 col, m_rows, m_cols);
    return m_pData[Ordering::IndexOf(row, col, m_rows, m_cols)];
  }

  const T &operator()(int row, int col) const {
    DEBUG_ASSERT(row < m_rows && col < m_cols,
                 "Out of bound array access ({} {}) with cols={}, rows={}", row,
                 col, m_rows, m_cols);
    return m_pData[Ordering::IndexOf(row, col, m_rows, m_cols)];
  }

  T *data() noexcept { return m_pData.data(); }

  const T *data() const noexcept { return m_pData.data(); }

  int rows() const noexcept { return m_rows; };
  int cols() const noexcept { return m_cols; };

  void reset() {
    for (size_t idx = 0; idx < m_pData.size(); ++idx) {
      m_pData[idx] = T();
    }
  }

  bool operator==(const Array &other) const {
    return m_rows == other.m_rows && m_cols == other.m_cols &&
           m_pData == other.m_pData;
  }

  bool operator!=(const Array &other) const { return !(*this == other); }

  /// Serialization

  void serialize(io::Serializer &serializer) const final {
    serializer.append(m_rows, "m_rows");
    serializer.append(m_cols, "m_cols");
    serializer.append(m_pData, "m_pData");
  }

  Array(io::Deserializer &deserializer)
      : m_rows(deserializer.restore<decltype(m_rows)>("m_rows")),
        m_cols(deserializer.restore<decltype(m_cols)>("m_cols")),
        m_pData(deserializer.restore<decltype(m_pData)>("m_pData")) {
    DEBUG_ASSERT(static_cast<int>(m_pData.size()) == m_rows * m_cols,
                 "Invalid size");
  }

private:
  int m_rows;
  int m_cols;
  std::vector<T> m_pData;
};
// Access is done in column major to enable possible vectorization
// by accessing the values of the same type contiguously, for example
// in memory we want XXXXXX.YYYYYY.ZZZZZZ....
using rArray = Array<double, RowMajor>;

/* The FixedArray class stores the values contiguously (all rows/cols in one
 * block).
 */
template <typename T, int _rows, int _cols, typename Ordering = RowMajor>
class FixedArray : public io::AbstractSerializable {
public:
  static const bool IsRowMajor = Ordering::IsRowMajor;

  explicit FixedArray(const T &defaultValue = T()) : m_pData{defaultValue} {}
  // support moving
  FixedArray(FixedArray &&rhs) = default;
  FixedArray &operator=(FixedArray &&rhs) = default;
  // support copying
  FixedArray(const FixedArray &rhs) = default;
  FixedArray &operator=(const FixedArray &rhs) = default;

  T &operator()(int row, int col) {
    DEBUG_ASSERT(row < _rows && col < _cols,
                 "Out of bound array access ({} {}) with cols={}, rows={}", row,
                 col, _rows, _cols);
    return m_pData[Ordering::IndexOf(row, col, _rows, _cols)];
  }

  const T &operator()(int row, int col) const {
    DEBUG_ASSERT(row < _rows && col < _cols,
                 "Out of bound array access ({} {}) with cols={}, rows={}", row,
                 col, _rows, _cols);
    return m_pData[Ordering::IndexOf(row, col, _rows, _cols)];
  }

  T *data() noexcept { return m_pData.data(); }

  const T *data() const noexcept { return m_pData.data(); }

  bool operator==(const FixedArray &other) const {
    return m_pData == other.m_pData;
  }

  bool operator!=(const FixedArray &other) const { return !(*this == other); }

  /// Serialization
  void serialize(io::Serializer &serializer) const final {
    serializer.append(m_pData, "m_pData");
  }

  FixedArray(io::Deserializer &deserializer)
      : m_pData(deserializer.restore<decltype(m_pData)>("m_pData")) {}

private:
  std::array<T, _rows * _cols> m_pData;
};
// Quick Matrix type
template <class T> using Matrix = FixedArray<T, 3, 3>;

} // namespace util

// TODO: use an actual linalg library like blaze/eigen/xtensor for this
template <class T> util::vec<T> operator-(const util::vec<T> &b) {
  return util::vec<T>(-b[0], -b[1], -b[2]);
}

namespace std {
std::ostream &operator<<(std::ostream &out, const util::rvec &v);
} // namespace std
#endif // ARRAY_H
