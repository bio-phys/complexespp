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
#ifndef ENERGY_PARAMETER_H
#define ENERGY_PARAMETER_H

#include "domains/beads.h"
#include "io/serializer.h"
#include "util/array.h"

namespace energy {

template <typename T> class PairParameter : public io::AbstractSerializable {
public:
  explicit PairParameter(util::Array<T> params) : m_params(std::move(params)) {}

  const T operator()(const domains::Bead a, const domains::Bead b) const
      noexcept {
    return m_params(a, b);
  }

  const T *data(const domains::Bead a) const {
    // It does not really matter if we are in
    // row major or column major as long as it is symmetric.
    return &m_params.data()[a * m_params.cols()];
  }

  size_t size() const {
    // Considering that rows = cols
    return m_params.rows();
  }

  bool operator==(const PairParameter<T> &other) const {
    return m_params == other.m_params;
  }

  bool operator!=(const PairParameter<T> &other) const {
    return !(*this == other);
  }

  void serialize(io::Serializer &serializer) const final {
    serializer.append(m_params, "m_params");
  }

  PairParameter(io::Deserializer &deserializer)
      : m_params(deserializer.restore<decltype(m_params)>("m_params")) {}

private:
  const util::Array<T> m_params;
};

} // namespace energy

#endif // ENERGY_PARAMETER_H
