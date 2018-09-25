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
#ifndef IO_FFPARAMETERS_H
#define IO_FFPARAMETERS_H

#include "energy/pairparameter.h"

namespace io {
energy::PairParameter<double>
readPairParameter(const std::string &file,
                  const std::vector<std::string> &beadTypes);

std::vector<std::string> readBeadTypes(const std::string &file);

std::vector<std::array<double, 8>>
readMembranePotential(const std::string &file,
                      const std::vector<std::string> &beadTypes);
} // namespace io

#endif // IO_FFPARAMETERS_H
