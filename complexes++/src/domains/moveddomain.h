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

  explicit MovedDomain(std::unique_ptr<AbstractDomain> &&inMovedDomain,
                       const std::string _type)
      : m_movedDomain(std::move(inMovedDomain)), m_type(_type) {}

  MovedDomain(const MovedDomain &) = delete;
  MovedDomain &operator=(const MovedDomain &) = delete;

public:
  MovedDomain() {}

  MovedDomain(MovedDomain &&) = default;
  MovedDomain &operator=(MovedDomain &&) = default;

  bool moveSucceed() const { return m_movedDomain != nullptr; }
  const std::string &type() const { return m_type; }

  std::unique_ptr<AbstractDomain> releaseDomain() {
    if (m_movedDomain == nullptr) {
      throw std::runtime_error(
          "Error MovedDomain::releaseDomain cannot be called on an "
          "unsuccessful MovedDomain");
    }
    return std::move(m_movedDomain);
  }

  static MovedDomain Fail() { return MovedDomain(); }

  static MovedDomain Success(std::unique_ptr<AbstractDomain> &&inMovedDomain,
                             const std::string &type) {
    return MovedDomain(std::move(inMovedDomain), type);
  }

  template <class DomainClass>
  static MovedDomain Success(DomainClass &&inMovedDomain,
                             const std::string &type) {
    return MovedDomain(
        std::unique_ptr<AbstractDomain>(
            std::move(std::make_unique<DomainClass>(std::move(inMovedDomain)))),
        type);
  }

  static MovedDomain Success(std::nullptr_t) = delete;
};
} // namespace domains

#endif // MOVEDDOMAIN_H
