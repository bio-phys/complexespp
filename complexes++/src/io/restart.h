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
#ifndef IO_RESTART_H
#define IO_RESTART_H

#include "util/random.h"

class string;

namespace mc {
class AbstractMcAlgo;
}

namespace pairkernels {
class PairKernelManager;
}

namespace energy {
class ForceField;
}

namespace io {
void saveRestart(const std::string &fname, const mc::AbstractMcAlgo &inAlgo,
                 const util::RNGEngine &inRng);
std::unique_ptr<mc::AbstractMcAlgo>
restoreRestart(const std::string &fname, const energy::ForceField &inForcefield,
               const pairkernels::PairKernelManager &inKernels,
               util::RNGEngine &inRng);
} // namespace io

#endif // IO_RESTART_H
