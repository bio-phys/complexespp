// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#ifndef IO_FFPARAMETERS_H
#define IO_FFPARAMETERS_H

#include "energy/pairparameter.h"

namespace io {
energy::PairParameter<double> readPairParameter(
    const std::string& file, const std::vector<std::string>& beadTypes);

std::vector<std::string> readBeadTypes(const std::string& file);

std::vector<std::array<double, 8>> readMembranePotential(
    const std::string& file, const std::vector<std::string>& beadTypes);
}  // namespace io

#endif  // IO_FFPARAMETERS_H
