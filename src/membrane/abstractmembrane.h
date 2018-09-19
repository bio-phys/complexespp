// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#ifndef MEMBRANE_ABSTRACTMEMBRANE_H
#define MEMBRANE_ABSTRACTMEMBRANE_H

#include "io/rebuilder.h"
#include "io/serializer.h"
#include "util/array.h"

namespace membrane {

template <typename T>
class AbstractMembrane : public io::RebuilderCore<AbstractMembrane<T>>,
                         public io::AbstractSerializable {
 public:
  virtual ~AbstractMembrane(){};
  virtual std::vector<T> distance(const util::rvec& bead,
                                  const util::rvec& box) const = 0;
  virtual util::Array<T> xyz() const = 0;
  virtual std::unique_ptr<AbstractMembrane<T>> copy() const = 0;
};

}  // namespace membrane
#endif  // MEMBRANE_ABSTRACTMEMBRANE_H
