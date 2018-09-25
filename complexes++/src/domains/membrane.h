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
#ifndef DOMAIN_MEMBRANEDOMAIN_H
#define DOMAIN_MEMBRANEDOMAIN_H

#include "domains/abstractdomain.h"
#include "domains/moveddomain.h"
#include "membrane/abstractmembrane.h"
#include "membrane/rebuildmembrane.h"

namespace domains {

class Membrane : public AbstractDomain {
public:
  explicit Membrane(
      const std::string &inTypename,
      std::unique_ptr<membrane::AbstractMembrane<double>> &&membrane,
      double z0_, double ePSI0_, int typeId_, int id_, Bead bead_,
      BeadChainID beadID)
      : AbstractDomain(inTypename, typeId_, id_, std::vector<Bead>{bead_},
                       std::vector<double>{0}, std::vector<BeadChainID>{beadID},
                       Connections(0)),
        m_membrane(std::move(membrane)), m_z0(z0_), m_ePSI0(ePSI0_) {
    setXyz(m_membrane->xyz());
  }

  MovedDomain move(const util::rvec &box, util::RNGEngine &rng) const final {
    UNUSED(box);
    UNUSED(rng);
    auto tmp = copy();
    return MovedDomain::Success(std::move(tmp), "membrane");
  }

  std::unique_ptr<AbstractDomain> copy() const final {
    auto c = Membrane(name(), m_membrane->copy(), m_z0, m_ePSI0, typeId(), id(),
                      beads()[0], BeadChainIDs()[0]);
    return std::make_unique<Membrane>(std::move(c));
  }

  std::vector<double> distance(const util::rvec &bead,
                               const util::rvec &box) const {
    return m_membrane->distance(bead, box);
  }

  static std::string Type() { return "membrane"; }

  virtual std::string type() const override { return Type(); }

  bool isMembrane() const final { return true; }
  double z0() const { return m_z0; }
  double ePSI0() const { return m_ePSI0; }

  // Serialization/Extraction

  void serialize(io::Serializer &serializer) const final {
    AbstractDomain::serializeCore(serializer);
    serializer.append(*m_membrane, "m_membrane");
    serializer.append(m_z0, "m_z0");
    serializer.append(m_ePSI0, "m_ePSI0");
  }

  Membrane(io::Deserializer &deserializer)
      : AbstractDomain(deserializer),
        m_membrane(
            membrane::RebuildMembrane<double>(deserializer, "m_membrane")),
        m_z0(deserializer.restore<decltype(m_z0)>("m_z0")),
        m_ePSI0(deserializer.restore<decltype(m_ePSI0)>("m_ePSI0")) {}

private:
  std::unique_ptr<membrane::AbstractMembrane<double>> m_membrane;
  double m_z0;
  double m_ePSI0;
};

} // namespace domains
#endif // DOMAIN_MEMBRANEDOMAIN_H
