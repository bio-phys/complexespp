// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#ifndef RIGIDDOMAIN_H
#define RIGIDDOMAIN_H

#include "domains/abstractdomain.h"
#include "domains/moveddomain.h"
#include "util/quaternions/quat.h"

namespace domains {

class Rigid : public AbstractDomain {
 public:
  explicit Rigid(const std::string& inTypename, int typeId_, int id_,
                 std::vector<Bead> beads_, std::vector<double> charges_,
                 std::vector<BeadChainID> beadChainIDs_,
                 const Connections& connections_, util::rvec trans, double phi);

  virtual MovedDomain move(const util::rvec& box,
                           util::RNGEngine& rng) const final;

  virtual std::unique_ptr<AbstractDomain> copy() const final;

  static std::string Type() { return "rigid"; }
  virtual std::string type() const final { return Type(); }

  void serialize(io::Serializer& serializer) const final {
    AbstractDomain::serializeCore(serializer);
    serializer.append(m_phi, "m_phi");
    serializer.append(m_trans, "m_trans");
  }

  Rigid(io::Deserializer& deserializer)
      : AbstractDomain(deserializer),
        m_phi(deserializer.restore<decltype(m_phi)>("m_phi")),
        m_trans(deserializer.restore<decltype(m_trans)>("m_trans")) {}

 private:
  double m_phi;
  util::rvec m_trans;
};
REBUILDER_REGISTER(Rigid);

}  // namespace domains

#endif  // RIGIDDOMAIN_H
