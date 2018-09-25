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
#ifndef IO_IO_H
#define IO_IO_H

#include "domains/system.h"
#include "io/reader.h"
#include "io/statfile.h"
#include "io/trajectoryfile.h"

namespace mc {
class AbstractMcAlgo;
}
class string;
class vector;
namespace setup {
class Config;
}
namespace energy {
class ForceField;
}
namespace pairkernels {
class PairKernelManager;
}

namespace io {

///////////////////////////
// Functions for Reading //
///////////////////////////
domains::System readTopologies(const std::string &yamlfile,
                               const std::vector<std::string> &beadTypes);
util::rvec readBox(const std::string &yamlfile);
energy::ForceField readForceField(const std::string &configPath,
                                  const setup::Config &config);
std::unique_ptr<mc::AbstractMcAlgo>
readRestart(const std::string &fname, const energy::ForceField &inForcefield,
            const pairkernels::PairKernelManager &inKernels,
            util::RNGEngine &inRng);
Reader readTrajectory(std::shared_ptr<domains::Domains> model,
                      const util::rvec &box, const std::string &file);
pairkernels::PairKernelManager
readKernelMapping(const std::string &yamlfile, const domains::Domains &domains);

///////////////////////////
// Functions for Writing //
///////////////////////////
void writeModel(TrajectoryFile &file, const domains::Domains &model,
                const util::rvec &box, const int i, const double time,
                const std::vector<std::string> &beadTypes);
void writeRestart(const std::string &fname, const mc::AbstractMcAlgo &inAlgo,
                  const util::RNGEngine inRng);

} // namespace io

#endif // IO_IO_H
