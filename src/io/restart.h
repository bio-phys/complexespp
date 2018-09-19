// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
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
void saveRestart(const std::string& fname, const mc::AbstractMcAlgo& inAlgo,
                 const util::RNGEngine& inRng);
std::unique_ptr<mc::AbstractMcAlgo> restoreRestart(
    const std::string& fname, const energy::ForceField& inForcefield,
    const pairkernels::PairKernelManager& inKernels, util::RNGEngine& inRng);
}  // namespace io

#endif  // IO_RESTART_H
