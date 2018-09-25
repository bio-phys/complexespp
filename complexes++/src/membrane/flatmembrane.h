// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#ifndef MEMBRANE_FLATMEMBRANE_H
#define MEMBRANE_FLATMEMBRANE_H

#include "membrane/abstractmembrane.h"
#include "util/pbc.h"

namespace membrane {

template <typename T>
class Flat : public AbstractMembrane<T> {
 public:
  explicit Flat(double zaxis) : m_zaxis(zaxis) {}

  std::vector<T> distance(const util::rvec& bead,
                          const util::rvec& box) const final {
    const auto diff = util::rvec(0, 0, bead[2] - m_zaxis);
    const auto dist = std::sqrt(util::pbc::DistSquare(diff, box));
    return std::vector<double>(1, dist);
  }

  util::Array<T> xyz() const final {
    auto pos = util::rArray(1, 3);
    pos(0, 2) = m_zaxis;
    return pos;
  }

  std::unique_ptr<AbstractMembrane<T>> copy() const final {
    auto tmp = *this;
    return std::make_unique<Flat<T>>(std::move(tmp));
  }
  static std::string Type() { return "flat"; }

  // Serialization/Extraction

  void serialize(io::Serializer& serializer) const final {
    serializer.append(Type(), "type");
    serializer.append(m_zaxis, "m_zaxis");
  }

  Flat(io::Deserializer& deserializer)
      : m_zaxis(deserializer.restore<decltype(m_zaxis)>("m_zaxis")) {}

 private:
  double m_zaxis;
};
using membrane_Flat_double = Flat<double>;
REBUILDER_REGISTER(membrane_Flat_double);

}  // namespace membrane
#endif  // MEMBRANE_FLATMEMBRANE_H
