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
#include <fmt/format.h>

#include "domains/beads.h"
#include "util/util.h"

namespace domains {

Bead findBeadID(const std::string &beadCode,
                const std::vector<std::string> &beadTypes) {
  const auto index =
      util::indexOf(std::cbegin(beadTypes), std::cend(beadTypes), beadCode);
  if (index == -1) {
    throw std::runtime_error(fmt::format("Unkown bead type: {}", beadCode));
  }

  return index;
}

} // namespace domains

namespace std {
std::ostream &operator<<(std::ostream &out, const domains::BeadChainID &bead) {
  out << fmt::format("{} {}", bead.beadID(), bead.chain());
  return out;
}
} // namespace std
