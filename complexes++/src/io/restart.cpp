// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------

#include <fstream>
#include <string>
#include <vector>

#include "io/restart.h"
#include "io/serializer.h"
#include "mc/abstractmcalgo.h"
#include "mc/simulation.h"
#include "util/random.h"
#include "util/util.h"
#include "version.h"

namespace io {

void saveRestart(const std::string& fname, const mc::AbstractMcAlgo& inAlgo,
                 const util::RNGEngine& inRng) {
#pragma omp critical(saveRestartFile)
  {
    auto serializer = Serializer();
    // write file header
    serializer.append(NAME, "program");
    serializer.append(GIT_COMMIT_HASH, "version");
    // write simulation state
    serializer.append(std::string(util::getRNGState(inRng)), "rng");
    serializer.append(inAlgo, "m_mcAlgo");

    auto buffer = serializer.getBuffer();
    std::ofstream out(fname.c_str(), std::ios::binary);
    out.write(reinterpret_cast<char*>(buffer.data()), buffer.size());
  }
}

std::vector<unsigned char> readFile(const std::string& filename) {
  std::ifstream file(filename.c_str(), std::ios::binary | std::ios::ate);
  const auto size = file.tellg();
  auto buffer = std::string(size, '\0');
  file.seekg(0);
  if (!file.read(&buffer[0], size)) {
    throw std::runtime_error("can't read restart file");
  }
  return std::vector<unsigned char>(std::begin(buffer), std::end(buffer));
}

std::unique_ptr<mc::AbstractMcAlgo> restoreRestart(
    const std::string& fname, const energy::ForceField& inForcefield,
    const pairkernels::PairKernelManager& inKernels, util::RNGEngine& inRng) {
  auto buffer = readFile(fname);
  auto deserializer = Deserializer(buffer.data(), buffer.size());

  // check file header for consistency
  const auto program = deserializer.restore<std::string>("program");
  if (program != NAME) {
    throw std::runtime_error("wrong restart file header\n");
  }
  const auto version = deserializer.restore<std::string>("version");
  if (version != GIT_COMMIT_HASH) {
    throw std::runtime_error("restart file from different version\n");
  }

  util::setRNGState(inRng, deserializer.restore<std::string>("rng"));
  auto algo = mc::Simulation::RebuildMcAlgo(deserializer, inForcefield,
                                            inKernels, inRng);
  algo->resetSessionCounter();
  return algo;
}

}  // namespace io
