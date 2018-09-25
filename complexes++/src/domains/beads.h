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
#ifndef DOMAINS_BEADS_H
#define DOMAINS_BEADS_H

#include <string>
#include <vector>

#include "io/serializer.h"
#include "util/util.h"

namespace domains {

using Bead = int;
constexpr auto NoneType = -1;

Bead findBeadID(const std::string &beadCode,
                const std::vector<std::string> &beadTypes);

class BeadChainID : public io::AbstractSerializable {
public:
  explicit BeadChainID(const std::string &chain, const Bead beadID)
      : m_chain(chain), m_beadID(beadID) {}

  const std::string chain() const { return m_chain; }
  int beadID() const { return m_beadID; }

  bool operator==(const BeadChainID &rhs) const {
    return m_chain == rhs.m_chain && m_beadID == rhs.m_beadID;
  }

  bool operator!=(const BeadChainID &rhs) const {
    return !(this->operator==(rhs));
  }

  void serialize(io::Serializer &serializer) const final {
    serializer.append(m_chain, "m_chain");
    serializer.append(m_beadID, "m_beadID");
  }

  BeadChainID(io::Deserializer &deserializer)
      : m_chain(deserializer.restore<decltype(m_chain)>("m_chain")),
        m_beadID(deserializer.restore<decltype(m_beadID)>("m_beadID")) {}

private:
  std::string m_chain;
  Bead m_beadID;
};

} // namespace domains

namespace std {
std::ostream &operator<<(std::ostream &out, const domains::BeadChainID &bead);
} // namespace std

#endif // DOMAINS_BEADS_H
