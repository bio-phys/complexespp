// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#ifndef IO_GRO_H
#define IO_GRO_H

#include "io/io.h"
#include <string>

namespace io {

void writeGRO(const std::string& out, const domains::Domains& model,
              const util::rvec& box, const std::vector<std::string>& beadTypes);

}  // namespace io

#endif  // IO_GRO_H
