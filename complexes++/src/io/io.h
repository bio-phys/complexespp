// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
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
domains::System readTopologies(const std::string& yamlfile,
                               const std::vector<std::string>& beadTypes);
util::rvec readBox(const std::string& yamlfile);
energy::ForceField readForceField(const std::string& configPath,
                                  const setup::Config& config);
std::unique_ptr<mc::AbstractMcAlgo> readRestart(
    const std::string& fname, const energy::ForceField& inForcefield,
    const pairkernels::PairKernelManager& inKernels, util::RNGEngine& inRng);
Reader readTrajectory(std::shared_ptr<domains::Domains> model,
                      const util::rvec& box, const std::string& file);
pairkernels::PairKernelManager readKernelMapping(
    const std::string& yamlfile, const domains::Domains& domains);

///////////////////////////
// Functions for Writing //
///////////////////////////
void writeModel(TrajectoryFile& file, const domains::Domains& model,
                const util::rvec& box, const int i, const double time,
                const std::vector<std::string>& beadTypes);
void writeRestart(const std::string& fname, const mc::AbstractMcAlgo& inAlgo,
                  const util::RNGEngine inRng);

}  // namespace io

#endif  // IO_IO_H
