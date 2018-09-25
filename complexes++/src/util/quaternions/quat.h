// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#ifndef UTIL_QUATERNIONS_QUAT_H
#define UTIL_QUATERNIONS_QUAT_H

#include "util/array.h"

namespace util {
namespace quaternions {

template <typename Real,
          typename = std::enable_if<std::is_floating_point<Real>::value>>
class Quaternion : public io::AbstractSerializable {
 public:
  Quaternion() : m_s(0), m_i(0), m_j(0), m_k(0) {}
  explicit Quaternion(const Real s, const Real i, const Real j, const Real k)
      : m_s(s), m_i(i), m_j(j), m_k(k) {
    DEBUG_ASSERT(std::abs(norm() - 1) < 1E-8,
                 "Quaternion is not normalized to 1, it is {}\n", norm());
  }
  explicit Quaternion(const util::vec<Real>& axis, const Real angle) {
    const auto r = std::sqrt(dot(axis, axis));
    m_s = cos(angle / 2);
    m_i = sin(angle / 2) * axis[0] / r;
    m_j = sin(angle / 2) * axis[1] / r;
    m_k = sin(angle / 2) * axis[2] / r;
  }
  // support moving
  Quaternion(Quaternion&& rhs) = default;
  Quaternion& operator=(Quaternion&& rhs) = default;
  // support copying
  Quaternion(const Quaternion& rhs) = default;
  Quaternion& operator=(const Quaternion& rhs) = default;

  double norm() const noexcept {
    return std::sqrt(m_s * m_s + m_i * m_i + m_j * m_j + m_k * m_k);
  }
  Quaternion<Real> conjugate() const noexcept {
    return Quaternion<Real>(m_s, -m_i, -m_j, -m_k);
  }
  util::vec<Real> vec() const noexcept {
    return util::vec<Real>(m_i, m_j, m_k);
  }

  Matrix<Real> toMat() const {
    auto mat = Matrix<Real>();
    // diagonal elements
    mat(0, 0) = m_s * m_s + m_i * m_i - m_j * m_j - m_k * m_k;
    mat(1, 1) = m_s * m_s - m_i * m_i + m_j * m_j - m_k * m_k;
    mat(2, 2) = m_s * m_s - m_i * m_i - m_j * m_j + m_k * m_k;
    // of diagonal elements
    mat(0, 1) = 2 * (m_i * m_j - m_s * m_k);
    mat(1, 0) = 2 * (m_i * m_j + m_s * m_k);

    mat(0, 2) = 2 * (m_i * m_k + m_s * m_j);
    mat(2, 0) = 2 * (m_i * m_k - m_s * m_j);

    mat(1, 2) = 2 * (m_j * m_k - m_s * m_i);
    mat(2, 1) = 2 * (m_j * m_k + m_s * m_i);

    return mat;
  }

  Quaternion<Real>& operator*=(const Quaternion<Real>& rhs) noexcept {
    const auto s =
        m_s * rhs.m_s - m_i * rhs.m_i - m_j * rhs.m_j - m_k * rhs.m_k;
    const auto i =
        m_s * rhs.m_i + m_i * rhs.m_s + m_j * rhs.m_k - m_k * rhs.m_j;
    const auto j =
        m_s * rhs.m_j - m_i * rhs.m_k + m_j * rhs.m_s + m_k * rhs.m_i;
    const auto k =
        m_s * rhs.m_k + m_i * rhs.m_j - m_j * rhs.m_i + m_k * rhs.m_s;
    m_s = s;
    m_i = i;
    m_j = j;
    m_k = k;
    return *this;
  }

  Quaternion<Real>& operator*=(const util::vec<Real>& v) noexcept {
    // when multiplying a vector and a Quaternion just treat the vector as a
    // Quaternion with a scalar value of 0.
    // NOTE: BE REALLY CAREFUL changing this!!
    const auto s = 0 - m_i * v[0] - m_j * v[1] - m_k * v[2];
    const auto i = m_s * v[0] + 0 + m_j * v[2] - m_k * v[1];
    const auto j = m_s * v[1] - m_i * v[2] + 0 + m_k * v[0];
    const auto k = m_s * v[2] + m_i * v[1] - m_j * v[0] + 0;
    m_s = s;
    m_i = i;
    m_j = j;
    m_k = k;
    return *this;
  }

  Quaternion<Real>& operator/=(const Real v) noexcept {
    m_s /= v;
    m_i /= v;
    m_j /= v;
    m_k /= v;
    return *this;
  }

  bool operator==(const Quaternion<Real>& other) const noexcept {
    return m_s == other.m_s || m_i == other.m_i || m_j == other.m_j ||
           m_k == other.m_k;
  }

  friend std::ostream& operator<<(std::ostream& out, const Quaternion& q) {
    auto constexpr format = "{:12.8f}{:12.8f}{:12.8f}{:12.8f}";
    out << fmt::format(format, q.m_s, q.m_i, q.m_j, q.m_k);
    return out;
  }

  std::array<Real, 4> rawValues() const { return {m_s, m_i, m_j, m_k}; }

  void serialize(io::Serializer& serializer) const final {
    serializer.append(m_s, "m_s");
    serializer.append(m_i, "m_i");
    serializer.append(m_j, "m_j");
    serializer.append(m_k, "m_k");
  }

  Quaternion(io::Deserializer& deserializer)
      : m_s(deserializer.restore<decltype(m_s)>("m_s")),
        m_i(deserializer.restore<decltype(m_i)>("m_i")),
        m_j(deserializer.restore<decltype(m_j)>("m_j")),
        m_k(deserializer.restore<decltype(m_s)>("m_k")) {}

 private:
  Real m_s;
  Real m_i;
  Real m_j;
  Real m_k;
};

template <typename Real,
          typename = std::enable_if<std::is_floating_point<Real>::value>>
constexpr Quaternion<Real> unity() noexcept {
  return Quaternion<Real>(1, 0, 0, 0);
}

template <typename Real,
          typename = std::enable_if<std::is_floating_point<Real>::value>>
Quaternion<Real> operator*(const Quaternion<Real>& lhs,
                           const Quaternion<Real>& rhs) {
  auto tmp = Quaternion<Real>(lhs);
  return tmp *= rhs;
}

template <typename Real,
          typename = std::enable_if<std::is_floating_point<Real>::value>>
Quaternion<Real> operator*(const Quaternion<Real>& lhs,
                           const vec<Real>& rhs) noexcept {
  auto tmp = Quaternion<Real>(lhs);
  return tmp *= rhs;
}

template <typename Real,
          typename = std::enable_if<std::is_floating_point<Real>::value>>
vec<Real> rotateVector(const Quaternion<Real> q, const vec<Real> v) noexcept {
  const auto tmp = q * v;
  const auto res = tmp * q.conjugate();
  return res.vec();
}

}  // namespace quaternions
}  // namespace util

#endif  // UTIL_QUATERNIONS_QUAT_H
