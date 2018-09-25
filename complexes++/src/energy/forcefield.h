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
#ifndef FORCEFIELD_H
#define FORCEFIELD_H

#include <string>
#include <vector>

#include "energy/pairparameter.h"
#include "util/file.h"

namespace energy {

class ForceField : public io::AbstractSerializable {
public:
  explicit ForceField(const std::vector<std::string> &_beadTypes,
                      const PairParameter<double> &_interActionEnergy,
                      const PairParameter<double> &_diameter,
                      const std::vector<double> &_chargeRadius,
                      const std::vector<std::array<double, 8>> &_membrane,
                      const double _debyeLength,
                      const double _dielectricConstant, const double _alpha)
      : m_beadTypes(_beadTypes), m_interActionEnergy(_interActionEnergy),
        m_diameter(_diameter), m_chargeRadius(_chargeRadius),
        m_membrane(_membrane), m_debyeLength(_debyeLength),
        m_dielectricConstant(_dielectricConstant), m_alpha(_alpha) {}

  const std::vector<std::string> &beadTypes() const { return m_beadTypes; }
  const PairParameter<double> &interActionEnergy() const {
    return m_interActionEnergy;
  }
  const PairParameter<double> &diameter() const { return m_diameter; }
  const std::vector<double> &chargeRadius() const { return m_chargeRadius; }
  const std::vector<std::array<double, 8>> &membrane() const {
    return m_membrane;
  }
  double debyeLength() const { return m_debyeLength; }
  double dielectricConstant() const { return m_dielectricConstant; }
  double alpha() const { return m_alpha; }

  //! Return true if the current object and other
  //! are equal (based on the comparison of the
  //! source directory)
  bool operator==(const ForceField &other) const {
    // clang-format off
    return m_debyeLength == other.m_debyeLength &&
           m_dielectricConstant == other.m_dielectricConstant &&
           m_alpha == other.m_alpha &&
           m_beadTypes == other.m_beadTypes &&
           m_interActionEnergy == other.m_interActionEnergy &&
           m_diameter == other.m_diameter &&
           m_chargeRadius == other.m_chargeRadius &&
           m_membrane == other.m_membrane;
    // clang-format on
  }

  bool operator!=(const ForceField &other) const { return !(*this == other); }

  //////////////////////////////////////////////////////////////////////////
  /// Serialize/deserialize
  //////////////////////////////////////////////////////////////////////////

  void serialize(io::Serializer &serializer) const final {
    serializer.append(m_beadTypes, "m_beadTypes");
    serializer.append(m_interActionEnergy, "m_interActionEnergy");
    serializer.append(m_diameter, "m_diameter");
    serializer.append(m_chargeRadius, "m_chargeRadius");
    serializer.append(m_membrane, "m_membrane");
    serializer.append(m_debyeLength, "m_debyeLength");
    serializer.append(m_dielectricConstant, "m_dielectricConstant");
    serializer.append(m_alpha, "m_alpha");
  }

  ForceField(io::Deserializer &deserializer)
      : m_beadTypes(deserializer.restore<decltype(m_beadTypes)>("m_beadTypes")),
        m_interActionEnergy(deserializer.restore<decltype(m_interActionEnergy)>(
            "m_interActionEnergy")),
        m_diameter(deserializer.restore<decltype(m_diameter)>("m_diameter")),
        m_chargeRadius(
            deserializer.restore<decltype(m_chargeRadius)>("m_chargeRadius")),
        m_membrane(deserializer.restore<decltype(m_membrane)>("m_membrane")),
        m_debyeLength(
            deserializer.restore<decltype(m_debyeLength)>("m_debyeLength")),
        m_dielectricConstant(
            deserializer.restore<decltype(m_dielectricConstant)>(
                "m_dielectricConstant")),
        m_alpha(deserializer.restore<decltype(m_alpha)>("m_alpha")) {}

private:
  /// Variable must remain constant to ensure equality for
  /// forcefield from the same directory
  std::vector<std::string> m_beadTypes;
  PairParameter<double> m_interActionEnergy;
  PairParameter<double> m_diameter;
  std::vector<double> m_chargeRadius;
  std::vector<std::array<double, 8>> m_membrane;
  double m_debyeLength;
  double m_dielectricConstant;
  double m_alpha;
};
} // namespace energy

#endif // FORCEFIELD_H
