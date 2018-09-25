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
#ifndef REBUILDMEMBRANE_HPP
#define REBUILDMEMBRANE_HPP

#include "membrane/abstractmembrane.h"
#include "membrane/flatmembrane.h"
#include "membrane/spheremembrane.h"
#include "membrane/tubemembrane.h"

namespace membrane {

template <typename T>
static std::unique_ptr<membrane::AbstractMembrane<T>>
RebuildMembrane(io::Deserializer &deserializer, const std::string key) {
  deserializer.access(key);
  return AbstractMembrane<T>::Rebuild(deserializer);
}
} // namespace membrane

#endif
