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

#include "domains/connection.h"
#include "constants.h"
#include "energy/energy.h"
#include "energy/forcefield.h"
#include "util/util.h"

namespace domains {

FlatConnection::FlatConnection() : Connection() {}
FlatConnection::FlatConnection(const int beadSelf_, const int domainId_,
                               const int beadOther_)
    : Connection(beadSelf_, domainId_, beadOther_) {}

double FlatConnection::energy(const double r2,
                              const energy::ForceField &forcefield) const {
  UNUSED(r2);
  UNUSED(forcefield);
  return 0;
}

GaussianConnection::GaussianConnection() : Connection(), m_N(1), m_k(.5) {}
GaussianConnection::GaussianConnection(const int beadSelf_, const int domainId_,
                                       const int beadOther_, const int N,
                                       const double b)
    : Connection(beadSelf_, domainId_, beadOther_), m_N(N),
      m_k(1.5 * (1 / b / b) / (N + 1)) {}

double GaussianConnection::energy(const double r2,
                                  const energy::ForceField &forcefield) const {
  UNUSED(forcefield);
  return m_k * r2;
}

HarmonicConnection::HarmonicConnection()
    : Connection(), m_x0(constants::polymerchain::bondLength),
      m_k(constants::polymerchain::kPseudoBond) {}
HarmonicConnection::HarmonicConnection(const int beadSelf_, const int domainId_,
                                       const int beadOther_, const double x0,
                                       const double k)
    : Connection(beadSelf_, domainId_, beadOther_), m_x0(x0), m_k(k) {}

double HarmonicConnection::energy(const double r2,
                                  const energy::ForceField &forcefield) const {
  UNUSED(forcefield);
  return energy::harmonic(std::sqrt(r2), m_x0, m_k);
}

} // namespace domains
