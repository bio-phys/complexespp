// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#ifndef MOVEDDOMAIN_H
#define MOVEDDOMAIN_H

#include <cstddef>
#include <memory>
#include <stdexcept>

namespace domains {
class AbstractDomain;

class MovedDomain {
  std::unique_ptr<AbstractDomain> m_movedDomain = nullptr;
  std::string m_type = "";

  explicit MovedDomain(std::unique_ptr<AbstractDomain>&& inMovedDomain,
                       const std::string _type)
      : m_movedDomain(std::move(inMovedDomain)), m_type(_type) {}

  MovedDomain(const MovedDomain&) = delete;
  MovedDomain& operator=(const MovedDomain&) = delete;

 public:
  MovedDomain() {}

  MovedDomain(MovedDomain&&) = default;
  MovedDomain& operator=(MovedDomain&&) = default;

  bool moveSucceed() const { return m_movedDomain != nullptr; }
  const std::string& type() const { return m_type; }

  std::unique_ptr<AbstractDomain> releaseDomain() {
    if (m_movedDomain == nullptr) {
      throw std::runtime_error(
          "Error MovedDomain::releaseDomain cannot be called on an "
          "unsuccessful MovedDomain");
    }
    return std::move(m_movedDomain);
  }

  static MovedDomain Fail() { return MovedDomain(); }

  static MovedDomain Success(std::unique_ptr<AbstractDomain>&& inMovedDomain,
                             const std::string& type) {
    return MovedDomain(std::move(inMovedDomain), type);
  }

  template <class DomainClass>
  static MovedDomain Success(DomainClass&& inMovedDomain,
                             const std::string& type) {
    return MovedDomain(
        std::unique_ptr<AbstractDomain>(
            std::move(std::make_unique<DomainClass>(std::move(inMovedDomain)))),
        type);
  }

  static MovedDomain Success(std::nullptr_t) = delete;
};
}  // namespace domains

#endif  // MOVEDDOMAIN_H
