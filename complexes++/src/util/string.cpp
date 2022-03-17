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
#include <string>
#include <vector>
#include <sstream>

#include "util/log.h"
#include "util/string.h"

namespace util {

std::vector<std::string> splitStr(const std::string &str,
                                  const char delimiter) {
	std::vector<std::string> tokens;
	std::string token;
	std::istringstream tokenStream(str);
	while (std::getline(tokenStream, token, delimiter)) {
		if (!token.empty()) {
			tokens.push_back(token);
		}
	}
	return tokens;
}

int truncateLeft(const int val, const int ndigits) {
  const int base = 10;
  const int max = static_cast<int>(std::pow(base, ndigits)) - 1;
  if (val <= max) {
    return val;
  }
  int mul = 1, leftDigit = val;
  while (leftDigit >= base) {
    leftDigit /= base;
    mul *= base;
  }
  return val - (leftDigit * mul);
}
} // namespace util
