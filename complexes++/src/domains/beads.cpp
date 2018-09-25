// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#include <fmt/format.h>

#include "domains/beads.h"
#include "util/util.h"

namespace domains {

Bead findBeadID(const std::string& beadCode,
                const std::vector<std::string>& beadTypes) {
  const auto index =
      util::indexOf(std::cbegin(beadTypes), std::cend(beadTypes), beadCode);
  if (index == -1) {
    throw std::runtime_error(fmt::format("Unkown bead type: {}", beadCode));
  }

  return index;
}

}  // namespace domains

namespace std {
std::ostream& operator<<(std::ostream& out, const domains::BeadChainID& bead) {
  out << fmt::format("{} {}", bead.beadID(), bead.chain());
  return out;
}
}  // namespace std
