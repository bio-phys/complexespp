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
#ifndef RIGIDDOMAIN_H
#define RIGIDDOMAIN_H

#include "domains/abstractdomain.h"
#include "domains/moveddomain.h"
#include "util/quaternions/quat.h"

namespace domains {

class Rigid : public AbstractDomain {
public:
  explicit Rigid(const std::string &inTypename, int typeId_, int id_,
                 std::vector<Bead> beads_, std::vector<double> charges_,
                 std::vector<BeadChainID> beadChainIDs_,
                 const Connections &connections_, util::rvec trans, double phi);

  virtual MovedDomain move(const util::rvec &box,
                           util::RNGEngine &rng) const final;

  virtual std::unique_ptr<AbstractDomain> copy() const final;

  static std::string Type() { return "rigid"; }
  virtual std::string type() const final { return Type(); }

  void serialize(io::Serializer &serializer) const final {
    AbstractDomain::serializeCore(serializer);
    serializer.append(m_phi, "m_phi");
    serializer.append(m_trans, "m_trans");
  }

  Rigid(io::Deserializer &deserializer)
      : AbstractDomain(deserializer),
        m_phi(deserializer.restore<decltype(m_phi)>("m_phi")),
        m_trans(deserializer.restore<decltype(m_trans)>("m_trans")) {}

private:
  double m_phi;
  util::rvec m_trans;
};
REBUILDER_REGISTER(Rigid);

} // namespace domains

#endif // RIGIDDOMAIN_H
