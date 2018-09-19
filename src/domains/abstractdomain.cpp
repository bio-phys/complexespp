// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------

#include <fmt/format.h>
#include <stdexcept>

#include "domains/abstractdomain.h"
#include "pairkernels/abstractpairkernel.h"
#include "util/pbc.h"

namespace domains {
// Passing by value here lets rvalues be moved and lvalue be copied. see EMC++
// Item 41
void AbstractDomain::setXyz(util::rArray xyz_) {
  if (xyz_.rows() != nBeads()) {
    throw std::invalid_argument(
        fmt::format("Number of beads in domain (id={}) does not match with "
                    "array dimensions: xyz.rows() = {}, expected {}\n",
                    id(), xyz_.rows(), nBeads()));
  }
  m_xyz = std::move(xyz_);
}

const util::rArray& AbstractDomain::xyz() const noexcept {
  return m_xyz;
}

util::rArray& AbstractDomain::xyz() noexcept {
  return m_xyz;
}

const std::vector<domains::Bead>& AbstractDomain::beads() const noexcept {
  return m_beads;
}

const std::vector<BeadChainID>& AbstractDomain::BeadChainIDs() const noexcept {
  return m_beadChainIDs;
}

const std::vector<double>& AbstractDomain::charges() const noexcept {
  return m_charges;
}

int AbstractDomain::id() const noexcept {
  return m_id;
}

const Connections& AbstractDomain::connections() const noexcept {
  return m_connections;
}

int AbstractDomain::nBeads() const noexcept {
  return static_cast<int>(m_beads.size());
}

const std::string& AbstractDomain::name() const noexcept {
  return m_name;
}

int AbstractDomain::typeId() const noexcept {
  return m_typeId;
}

void AbstractDomain::energyForAllConnections(
    const domains::Domains& others, const util::rvec& box,
    const energy::ForceField& forcefield, energy::rEnergyMatrix& outRes) const {
  const auto& xyzThis = xyz();
  auto d = util::rvec();

  for (const auto& con : m_connections) {
    const auto i = con->beadSelf();
    const auto j = con->beadOther();

    const auto idOther = con->domainId();
    const auto& xyzOther = others[idOther]->xyz();

    d[0] = xyzThis(i, 0) - xyzOther(j, 0);
    d[1] = xyzThis(i, 1) - xyzOther(j, 1);
    d[2] = xyzThis(i, 2) - xyzOther(j, 2);
    const double energyRes =
        con->energy(util::pbc::DistSquare(d, box), forcefield);
    outRes.addEnergyConnections(0, idOther, energyRes);
  }
}

double AbstractDomain::energyForConnections(
    const AbstractDomain& other, const util::rvec& box,
    const energy::ForceField& forcefield) const {
  if (id() == other.id()) {
    return 0;
  }

  const auto& xyzOther = other.xyz();
  const auto& xyzThis = xyz();
  auto d = util::rvec();

  // The domains are sorted, to find the interval corresponding to other.id()
  // we can use a binary search of the extremities (bounds)

  const auto firstCon = std::lower_bound(
      m_connections.begin(), m_connections.end(), other.id(),
      [](const std::shared_ptr<Connection>& con, const int idOther) {
        return con->domainId() < idOther;
      });

  const auto lastCon = std::upper_bound(
      firstCon, m_connections.end(), other.id(),
      [](const int idOther, const std::shared_ptr<Connection>& con) {
        return idOther < con->domainId();
      });

  // Loop over connections
  auto enLink = 0.0;
  for (auto iter = firstCon; iter != lastCon; ++iter) {
    const auto& con = (*iter);
    // Could use assert(other.id() == con->domainId());
    const auto i = con->beadSelf();
    const auto j = con->beadOther();
    d[0] = xyzThis(i, 0) - xyzOther(j, 0);
    d[1] = xyzThis(i, 1) - xyzOther(j, 1);
    d[2] = xyzThis(i, 2) - xyzOther(j, 2);
    enLink += con->energy(util::pbc::DistSquare(d, box), forcefield);
  }
  return enLink;
}

template <typename T>
void checkSize(int n, std::vector<T> col, const std::string& name, int id) {
  if (col.size() != static_cast<std::size_t>(n)) {
    throw std::invalid_argument(
        fmt::format("Number of beads in domain (id={}) does not match with "
                    "vector size: {}.size() = {}, expected {}\n",
                    id, name, col.size(), n));
  }
}

AbstractDomain::AbstractDomain(const std::string& inName, int typeId_, int id_,
                               std::vector<Bead> beads_,
                               std::vector<double> charges_,
                               std::vector<BeadChainID> beadChainIDs_,
                               const Connections& connections_)
    : m_id(id_),
      m_name(inName),
      m_typeId(typeId_),
      m_xyz(static_cast<int>(beads_.size()), 3),
      m_beads(beads_),
      m_charges(charges_),
      m_beadChainIDs(beadChainIDs_),
      m_connections(connections_) {
  checkSize(nBeads(), m_charges, "charges", id());
  checkSize(nBeads(), m_beadChainIDs, "beadChainIDs", id());

  // Sort the domains based on others' ids
  std::sort(m_connections.begin(), m_connections.end(),
            [](const std::shared_ptr<Connection>& c1,
               const std::shared_ptr<Connection>& c2) {
              return c1->domainId() < c2->domainId();
            });
}
}  // namespace domains
