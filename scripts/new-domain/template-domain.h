// clang-format off
#ifndef @UPPERCASEDOMAIN@_H
#define @UPPERCASEDOMAIN@_H

//! This domain has been generated automatically
//! by @USER@ the @DATE@

#include "domains/abstractdomain.h"

namespace domains {

class @CLASSNAME@ : public Domain {
 public:
  explicit @CLASSNAME@(const std::string& inTypename, int typeId_, int id_,
                 std::vector<Bead> beads_, std::vector<double> charges_,
                 std::vector<BeadChainID> beadChainIDs_,
                 const Connections& connections_)
    : Domain(inTypename, typeId_, id_, beads_, charges_, beadChainIDs_,
             connections_) {}

  double energyInner(const util::rvec& box,
                     const energy::ForceField& forcefield) const final {
    // Notice: decide how to compute inner energy
    static_assert(false,
                  "Implement this method or remove this assert"
                  "if nothing has to be done here.");
    return 0;
  }

  MovedDomain move(const util::rvec& box, const Domains& topology,
                   util::RNGEngine& rng,
                   const AbstractCollisionDetector<double>& collisionDetector) const final {
    // Notice: decide how a domain should move
  }

  static std::string Type() { return "@STRINGDOMAIN@"; }
  std::string type() const final {
    return Type();
  }

  std::unique_ptr<Domain> copy() const final {
    auto tmp = *this;
    return std::make_unique<@CLASSNAME@>(std::move(tmp));
  }

  void serialize(io::Serializer& serializer) const final {
    Domain::serializeCore(serializer);
  }

  #CLASSNAME@(io::Deserializer& deserializer)
    : Domain(deserializer) {}

 private:
  // Notice: Put the class attributes here
};
REBUILDER_REGISTER(@CLASSNAME@);
} // namespace domains

#endif // @UPPERCASEDOMAIN@_H
// clang-format on
