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
#ifndef UTIL_STRING_H
#define UTIL_STRING_H

class string;
class vector;

namespace util {

std::vector<std::string> splitStr(const std::string &str,
                                  const char delimiter);

/* truncate left most number from an integer. So 12345 -> 2345
 * Useful to prevent overflow in writing a PDB/GRO
 */
int truncateLeft(const int val, const int ndigits = 5);

} // namespace util

#endif // UTIL_STRING_H
