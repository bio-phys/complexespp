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
#ifndef ABSTRACTDOMAIN_H
#define ABSTRACTDOMAIN_H

#include <memory>
#include <utility>
#include <vector>

#include "domains/beads.h"
#include "domains/connection.h"
#include "domains/moveddomain.h"
#include "energy/energymatrix.h"
#include "io/rebuilder.h"
#include "io/serializer.h"
#include "util/array.h"
#include "util/random.h"

namespace energy {
class Forcefield;
}

namespace energy {
class AbstractPairKernel;
class ComputeResult;
} // namespace energy

namespace domains {

// make type aliases available at the header beginning.
class AbstractDomain;
using Domains = std::vector<std::unique_ptr<AbstractDomain>>;

class AbstractDomain : public io::RebuilderCore<AbstractDomain>,
                       public io::AbstractSerializable {
public:
  AbstractDomain(const std::string &inTypename, int typeId_, int id_,
                 std::vector<Bead> beads_, std::vector<double> charges_,
                 std::vector<BeadChainID> beadChainIDs_,
                 const Connections &connections_);
  virtual ~AbstractDomain() {}

  // support moves
  explicit AbstractDomain(AbstractDomain &&rhs) = default;
  AbstractDomain &operator=(AbstractDomain &&rhs) = default;
  // support copying
  AbstractDomain(const AbstractDomain &rhs) = default;
  AbstractDomain &operator=(const AbstractDomain &rhs) = default;

  /** This function computes the energy for the connections
   * between the current AbstractDomain and the other AbstractDomain.
   * This function is called once. If it is called then
   * energyForAllConnections will certainly not be called.
   */
  double energyForConnections(const AbstractDomain &other,
                              const util::rvec &box,
                              const energy::ForceField &forcefield) const;

  /** This function computes the energy for the connections
   * between the current AbstractDomain and all the other domains.
   * This function is called once. If it is called then
   * energyForConnections will certainly not be called. AND SHOULDN'T !!
   * The expected result is equal to:
   * @code for(all domain d in others)
   * @code    outRes += energyForConnections(d, box, forcefield)
   * @code endfor
   */
  void energyForAllConnections(const domains::Domains &others,
                               const util::rvec &box,
                               const energy::ForceField &forcefield,
                               energy::rEnergyMatrix &outRes) const;
  // API to define by child
  virtual MovedDomain move(const util::rvec &box,
                           util::RNGEngine &rng) const = 0;
  virtual std::unique_ptr<AbstractDomain> copy() const = 0;
  virtual std::string type() const = 0;
  // Membrane stuff (refactor out of this when possible)
  virtual bool isMembrane() const { return false; }

  // Make getting the values easy.
  const util::rArray &xyz() const noexcept;
  util::rArray &xyz() noexcept;
  const std::vector<domains::Bead> &beads() const noexcept;
  const std::vector<BeadChainID> &BeadChainIDs() const noexcept;
  const std::vector<double> &charges() const noexcept;
  int id() const noexcept;
  const Connections &connections() const noexcept;
  const std::string &name() const noexcept;
  int typeId() const noexcept;
  int nBeads() const noexcept;

  // Setter
  void setXyz(util::rArray xyz_);

  static Connections RebuildConnections(io::Deserializer &deserializer) {
    Connections connections;
    size_t m_connections_size =
        deserializer.restore<decltype(connections.size())>(
            "m_connections.size()");
    connections.resize(m_connections_size);
    for (size_t idx_connections = 0; idx_connections < m_connections_size;
         ++idx_connections) {
      deserializer.access("connection");
      connections[idx_connections] = Connection::Rebuild(deserializer);
    }

    return connections;
  }

  AbstractDomain(io::Deserializer &deserializer)
      : m_id(deserializer.restore<decltype(m_id)>("m_id")),
        m_name(deserializer.restore<decltype(m_name)>("m_name")),
        m_typeId(deserializer.restore<decltype(m_typeId)>("m_typeId")),
        m_xyz(deserializer.restore<decltype(m_xyz)>("m_xyz")),
        m_beads(deserializer.restore<decltype(m_beads)>("m_beads")),
        m_charges(deserializer.restore<decltype(m_charges)>("m_charges")),
        m_beadChainIDs(
            deserializer.restore<decltype(m_beadChainIDs)>("m_beadChainIDs")),
        m_connections(RebuildConnections(deserializer)) {}

protected:
  void serializeCore(io::Serializer &serializer) const {
    serializer.append(type(), "type");
    serializer.append(m_id, "m_id");
    serializer.append(m_name, "m_name");
    serializer.append(m_typeId, "m_typeId");
    serializer.append(m_xyz, "m_xyz");
    serializer.append(m_beads, "m_beads");
    serializer.append(m_charges, "m_charges");
    serializer.append(m_beadChainIDs, "m_beadChainIDs");
    serializer.append(m_connections.size(), "m_connections.size()");
    for (size_t idx_connections = 0; idx_connections < m_connections.size();
         ++idx_connections) {
      serializer.append(*m_connections[idx_connections], "connection");
    }
  }

private:
  int m_id;
  std::string m_name;
  int m_typeId;
  util::rArray m_xyz;
  std::vector<Bead> m_beads;
  std::vector<double> m_charges;
  std::vector<BeadChainID> m_beadChainIDs;
  Connections m_connections;
};
} // namespace domains

#endif // ABSTRACTDOMAIN_H
