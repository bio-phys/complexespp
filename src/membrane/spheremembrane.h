// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#ifndef MEMBRANE_SPHEREMEMBRANE_H
#define MEMBRANE_SPHEREMEMBRANE_H

#include "membrane/abstractmembrane.h"
#include "util/pbc.h"

namespace membrane {

template <typename T>
class Sphere : public AbstractMembrane<T> {
 public:
  explicit Sphere(double radius, util::rvec center)
      : m_radius(radius), m_center(center) {}

  std::vector<T> distance(const util::rvec& bead,
                          const util::rvec& box) const final {
    const auto diff = util::rvec(bead[0] - m_center[0], bead[1] - m_center[1],
                                 bead[2] - m_center[2]);
    const auto dist =
        std::abs(std::sqrt(util::pbc::DistSquare(diff, box)) - m_radius);
    return std::vector<double>(1, dist);
  }

  util::Array<T> xyz() const final {
    auto pos = util::rArray(1, 3);
    pos(0, 0) = m_center[0];
    pos(0, 1) = m_center[1];
    pos(0, 2) = m_center[2];
    return pos;
  }

  std::unique_ptr<AbstractMembrane<T>> copy() const final {
    auto tmp = *this;
    return std::make_unique<Sphere<T>>(std::move(tmp));
  }
  static std::string Type() { return "sphere"; }

  // Serialization/Extraction

  void serialize(io::Serializer& serializer) const final {
    serializer.append(Type(), "type");
    serializer.append(m_radius, "m_radius");
    serializer.append(m_center, "m_center");
  }

  Sphere(io::Deserializer& deserializer)
      : m_radius(deserializer.restore<decltype(m_radius)>("m_radius")),
        m_center(deserializer.restore<decltype(m_center)>("m_center")) {}

 private:
  double m_radius;
  util::rvec m_center;
};
using membrane_Sphere_double = Sphere<double>;
REBUILDER_REGISTER(membrane_Sphere_double);

}  // namespace membrane
#endif  // MEMBRANE_SPHEREMEMBRANE_H
