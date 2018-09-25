// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#include <fmt/format.h>

#include "constants.h"
#include "domains/abstractdomain.h"
#include "io/cplx.h"
#include "io/ffparameters.h"
#include "io/gro.h"
#include "io/io.h"
#include "io/pdb.h"
#include "io/restart.h"
#include "io/xdr.h"
#include "mc/abstractmcalgo.h"
#include "setup/config.h"
#include "util/file.h"
#include "util/log.h"
#include "version.h"

namespace io {

domains::System readTopologies(const std::string& file,
                               const std::vector<std::string>& beadTypes) {
  util::throwIfFileDoesNotExists(file);
  const auto suffix = util::fileSuffix(file);

  if (suffix == ".cplx") {
    return readCPLX(file, beadTypes);
  } else {
    throw std::invalid_argument(
        fmt::format("{} <-- unkown file format.\n", file));
  }
}

util::rvec readBox(const std::string& file) {
  util::throwIfFileDoesNotExists(file);
  const auto suffix = util::fileSuffix(file);

  if (suffix == ".cplx") {
    return readBoxCPLX(file);
  } else {
    throw std::invalid_argument(
        fmt::format("{} <-- unkown file format.\n", file));
  }
}

energy::ForceField readForceField(const std::string& configPath,
                                  const setup::Config& config) {
  return readForceFieldCPLX(util::appendPathsIfSubRelative(
      configPath, config.value<std::string>("structure")));
}

pairkernels::PairKernelManager readKernelMapping(
    const std::string& file, const domains::Domains& domains) {
  util::throwIfFileDoesNotExists(file);
  const auto suffix = util::fileSuffix(file);

  if (suffix == ".cplx") {
    return readCPLXKernels(file, domains);
  } else {
    throw std::invalid_argument(
        fmt::format("{} <-- unkown file format.\n", file));
  }
}

std::unique_ptr<mc::AbstractMcAlgo> readRestart(
    const std::string& file, const energy::ForceField& inForcefield,
    const pairkernels::PairKernelManager& inKernels, util::RNGEngine& inRng) {
  util::throwIfFileDoesNotExists(file);
  return restoreRestart(file, inForcefield, inKernels, inRng);
}

void writeModel(TrajectoryFile& file, const domains::Domains& model,
                const util::rvec& box, const int i, const double time,
                const std::vector<std::string>& beadTypes) {
  if (file.type() == "pdb") {
    writePDB(file.fstream(), model, box, i, beadTypes);
  } else if (file.type() == "xtc") {
    writeXTC(file.xdr(), model, box, i, time);
  } else if (file.type() == "trr") {
    writeTRR(file.xdr(), model, box, i, time);
  } else if (file.type() == "gro") {
    writeGRO(file.fname(), model, box, beadTypes);
  } else {
    throw std::invalid_argument(
        fmt::format("{} <-- unkown file format.\n", file.fname()));
  }
}

void writeRestart(const std::string& fname, const mc::AbstractMcAlgo& inAlgo,
                  const util::RNGEngine inRng) {
  saveRestart(fname, inAlgo, inRng);
}

Reader readTrajectory(std::shared_ptr<domains::Domains> dom,
                      const util::rvec& box, const std::string& file) {
  util::throwIfFileDoesNotExists(file);
  const auto suffix = util::fileSuffix(file);

  if (suffix == ".pdb") {
    return std::make_unique<PDBReader>(dom, box, file);
  } else if (suffix == ".xtc" || suffix == ".trr") {
    return std::make_unique<XDRReader>(dom, box, file);
  } else {
    throw std::invalid_argument(
        fmt::format("{} <-- unkown file format.\n", file));
  }
}
}  // namespace io
