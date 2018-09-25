// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#ifndef DOMAINS_BEADS_H
#define DOMAINS_BEADS_H

#include <string>
#include <vector>

#include "io/serializer.h"
#include "util/util.h"

namespace domains {

using Bead = int;
constexpr auto NoneType = -1;

Bead findBeadID(const std::string& beadCode,
                const std::vector<std::string>& beadTypes);

class BeadChainID : public io::AbstractSerializable {
 public:
  explicit BeadChainID(const std::string& chain, const Bead beadID)
      : m_chain(chain), m_beadID(beadID) {}

  const std::string chain() const { return m_chain; }
  int beadID() const { return m_beadID; }

  bool operator==(const BeadChainID& rhs) const {
    return m_chain == rhs.m_chain && m_beadID == rhs.m_beadID;
  }

  bool operator!=(const BeadChainID& rhs) const {
    return !(this->operator==(rhs));
  }

  void serialize(io::Serializer& serializer) const final {
    serializer.append(m_chain, "m_chain");
    serializer.append(m_beadID, "m_beadID");
  }

  BeadChainID(io::Deserializer& deserializer)
      : m_chain(deserializer.restore<decltype(m_chain)>("m_chain")),
        m_beadID(deserializer.restore<decltype(m_beadID)>("m_beadID")) {}

 private:
  std::string m_chain;
  Bead m_beadID;
};

}  // namespace domains

namespace std {
std::ostream& operator<<(std::ostream& out, const domains::BeadChainID& bead);
}  // namespace std

#endif  // DOMAINS_BEADS_H
