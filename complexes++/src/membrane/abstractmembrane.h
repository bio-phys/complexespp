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
  virtual std::vector<T> distance(const util::rvec &bead,
                                  const util::rvec &box) const = 0;
  virtual util::Array<T> xyz() const = 0;
  virtual std::unique_ptr<AbstractMembrane<T>> copy() const = 0;
};

} // namespace membrane
#endif // MEMBRANE_ABSTRACTMEMBRANE_H
