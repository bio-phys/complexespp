// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#ifndef DOMAINS_CONNECTION_H
#define DOMAINS_CONNECTION_H

#include <memory>
#include <vector>

#include "io/rebuilder.h"
#include "io/serializer.h"

namespace energy {
class ForceField;
}  // namespace energy

namespace domains {
// make type aliases available at the header beginning.
class Connection;
using Connections = std::vector<std::shared_ptr<Connection>>;

class Connection : public io::RebuilderCore<Connection>,
                   public io::AbstractSerializable {
 public:
  explicit Connection() : m_beadSelf(-1), m_beadOther(-1), m_domainId(-1) {}
  explicit Connection(const int beadSelf, const int domainId,
                      const int beadOther)
      : m_beadSelf(beadSelf), m_beadOther(beadOther), m_domainId(domainId) {}
  virtual ~Connection() {}

  bool operator==(const Connection& rhs) const {
    return m_beadSelf == rhs.m_beadSelf && m_beadOther == rhs.m_beadOther &&
           m_domainId == rhs.m_domainId;
  }

  bool operator!=(const Connection& rhs) const {
    return !(this->operator==(rhs));
  }

  virtual double energy(const double r2,
                        const energy::ForceField& forcefield) const = 0;

  int beadSelf() const { return m_beadSelf; }
  int beadOther() const { return m_beadOther; }
  int domainId() const { return m_domainId; }

  virtual std::string type() const = 0;

  Connection(io::Deserializer& deserializer)
      : m_beadSelf(deserializer.restore<decltype(m_beadSelf)>("m_beadSelf")),
        m_beadOther(deserializer.restore<decltype(m_beadOther)>("m_beadOther")),
        m_domainId(deserializer.restore<decltype(m_domainId)>("m_domainId")) {}

 protected:
  virtual void serializeCore(io::Serializer& serializer) const {
    serializer.append(type(), "type");
    serializer.append(m_beadSelf, "m_beadSelf");
    serializer.append(m_beadOther, "m_beadOther");
    serializer.append(m_domainId, "m_domainId");
  }

 private:
  int m_beadSelf;
  int m_beadOther;
  int m_domainId;
};

class FlatConnection : public Connection {
 public:
  explicit FlatConnection();
  explicit FlatConnection(const int beadSelf_, const int domainId_,
                          const int beadOther_);
  double energy(const double r2,
                const energy::ForceField& forcefield) const final;

  static std::string Type() { return "FlatConnection"; }
  std::string type() const final { return Type(); }

  void serialize(io::Serializer& serializer) const final {
    Connection::serializeCore(serializer);
  }
  FlatConnection(io::Deserializer& deserializer) : Connection(deserializer) {}
};
REBUILDER_REGISTER(FlatConnection);

class GaussianConnection : public Connection {
 public:
  explicit GaussianConnection();
  explicit GaussianConnection(const int beadSelf_, const int domainId_,
                              const int beadOther_, const int N,
                              const double b);
  double energy(const double r2,
                const energy::ForceField& forcefield) const final;

  static std::string Type() { return "GaussianConnection"; }
  std::string type() const final { return Type(); }

  void serialize(io::Serializer& serializer) const final {
    Connection::serializeCore(serializer);
    serializer.append(m_N, "m_N");
    serializer.append(m_k, "m_k");
  }

  GaussianConnection(io::Deserializer& deserializer)
      : Connection(deserializer),
        m_N(deserializer.restore<decltype(m_N)>("m_N")),
        m_k(deserializer.restore<decltype(m_k)>("m_k")) {}

 private:
  int m_N;
  double m_k;
};
REBUILDER_REGISTER(GaussianConnection);

class HarmonicConnection : public Connection {
 public:
  explicit HarmonicConnection();
  explicit HarmonicConnection(const int beadSelf_, const int domainId_,
                              const int beadOther_, const double x0,
                              const double k);
  double energy(const double r2,
                const energy::ForceField& forcefield) const final;

  static std::string Type() { return "HarmonicConnection"; }
  std::string type() const final { return Type(); }

  void serialize(io::Serializer& serializer) const final {
    Connection::serializeCore(serializer);
    serializer.append(m_x0, "m_x0");
    serializer.append(m_k, "m_k");
  }

  HarmonicConnection(io::Deserializer& deserializer)
      : Connection(deserializer),
        m_x0(deserializer.restore<decltype(m_x0)>("m_x0")),
        m_k(deserializer.restore<decltype(m_k)>("m_k")) {}

 private:
  double m_x0;
  double m_k;
};
REBUILDER_REGISTER(HarmonicConnection);

}  // namespace domains

#endif  // DOMAINS_CONNECTION_H
