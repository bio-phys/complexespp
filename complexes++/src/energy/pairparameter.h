// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#ifndef ENERGY_PARAMETER_H
#define ENERGY_PARAMETER_H

#include "domains/beads.h"
#include "io/serializer.h"
#include "util/array.h"

namespace energy {

template <typename T>
class PairParameter : public io::AbstractSerializable {
 public:
  explicit PairParameter(util::Array<T> params) : m_params(std::move(params)) {}

  const T operator()(const domains::Bead a, const domains::Bead b) const
      noexcept {
    return m_params(a, b);
  }

  const T* data(const domains::Bead a) const {
    // It does not really matter if we are in
    // row major or column major as long as it is symmetric.
    return &m_params.data()[a * m_params.cols()];
  }

  size_t size() const {
    // Considering that rows = cols
    return m_params.rows();
  }

  bool operator==(const PairParameter<T>& other) const {
    return m_params == other.m_params;
  }

  bool operator!=(const PairParameter<T>& other) const {
    return !(*this == other);
  }

  void serialize(io::Serializer& serializer) const final {
    serializer.append(m_params, "m_params");
  }

  PairParameter(io::Deserializer& deserializer)
      : m_params(deserializer.restore<decltype(m_params)>("m_params")) {}

 private:
  const util::Array<T> m_params;
};

}  // namespace energy

#endif  // ENERGY_PARAMETER_H
