// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#ifndef MEMBRANE_TUBEMEMBRANE_H
#define MEMBRANE_TUBEMEMBRANE_H

#include "membrane/abstractmembrane.h"
#include "util/pbc.h"

namespace membrane {

template <typename T>
class Tube : public AbstractMembrane<T> {
 public:
  explicit Tube(double x, double y, double radius)
      : m_radius(radius), m_x(x), m_y(y) {}

  // Only works when whole tube is inside of box and doesn't distinguish between
  // inside and outside of tube.
  std::vector<T> distance(const util::rvec& bead,
                          const util::rvec& box) const final {
    const auto diff = util::rvec(bead[0] - m_x, bead[1] - m_y, 0);
    const auto dist =
        std::abs(std::sqrt(util::pbc::DistSquare(diff, box)) - m_radius);
    return std::vector<double>(1, dist);
  }

  util::Array<T> xyz() const final {
    auto pos = util::rArray(1, 3);
    pos(0, 0) = m_x;
    pos(0, 1) = m_y;
    return pos;
  }

  std::unique_ptr<AbstractMembrane<T>> copy() const final {
    auto tmp = *this;
    return std::make_unique<Tube>(std::move(tmp));
  }
  static std::string Type() { return "tube"; }

  // Serialization/Extraction

  void serialize(io::Serializer& serializer) const final {
    serializer.append(Type(), "type");
    serializer.append(m_radius, "m_radius");
    serializer.append(m_x, "m_x");
    serializer.append(m_y, "m_y");
  }

  Tube(io::Deserializer& deserializer)
      : m_radius(deserializer.restore<decltype(m_radius)>("m_radius")),
        m_x(deserializer.restore<decltype(m_x)>("m_x")),
        m_y(deserializer.restore<decltype(m_y)>("m_y")) {}

 private:
  double m_radius;
  double m_x;
  double m_y;
};
using membrane_Tube_double = Tube<double>;
REBUILDER_REGISTER(membrane_Tube_double);

}  // namespace membrane
#endif  // MEMBRANE_TUBEMEMBRANE_H
