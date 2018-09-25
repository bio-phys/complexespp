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
#ifndef ABSTRACTEXCHANGEBUILDER_H
#define ABSTRACTEXCHANGEBUILDER_H

#include <string>

namespace mc {

//! This abstract class represents the interface
//! to a high level exchange alorithm.
//! The current usage is to call:
//! - init : to load from disk and restart if needed
//! - compareSimulations : to ensure compatible configs
//! - run : to execute the complete exchange algo
class AbstractExchangeSimulation {
public:
  virtual ~AbstractExchangeSimulation() {}
  virtual bool compareSimulations(const bool assertOnFailure) = 0;
  virtual int run() = 0;
};
} // namespace mc

#endif
