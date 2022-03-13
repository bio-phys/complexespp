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
#include <array>
#include <cstdint>
#include <iostream>
#include <vector>

#include "util/array.h"

namespace std {
std::ostream &operator<<(std::ostream &out, const util::rvec &v) {
  out << fmt::format("{:12.8f},{:12.8f},{:12.8f}", v[0], v[1], v[2]);
  return out;
} 
} // namespace std
