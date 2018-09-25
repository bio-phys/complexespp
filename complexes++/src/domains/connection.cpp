// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------

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
                              const energy::ForceField& forcefield) const {
  UNUSED(r2);
  UNUSED(forcefield);
  return 0;
}

GaussianConnection::GaussianConnection() : Connection(), m_N(1), m_k(.5) {}
GaussianConnection::GaussianConnection(const int beadSelf_, const int domainId_,
                                       const int beadOther_, const int N,
                                       const double b)
    : Connection(beadSelf_, domainId_, beadOther_),
      m_N(N),
      m_k(1.5 * (1 / b / b) / (N + 1)) {}

double GaussianConnection::energy(const double r2,
                                  const energy::ForceField& forcefield) const {
  UNUSED(forcefield);
  return m_k * r2;
}

HarmonicConnection::HarmonicConnection()
    : Connection(),
      m_x0(constants::polymerchain::bondLength),
      m_k(constants::polymerchain::kPseudoBond) {}
HarmonicConnection::HarmonicConnection(const int beadSelf_, const int domainId_,
                                       const int beadOther_, const double x0,
                                       const double k)
    : Connection(beadSelf_, domainId_, beadOther_), m_x0(x0), m_k(k) {}

double HarmonicConnection::energy(const double r2,
                                  const energy::ForceField& forcefield) const {
  UNUSED(forcefield);
  return energy::harmonic(std::sqrt(r2), m_x0, m_k);
}

}  // namespace domains
