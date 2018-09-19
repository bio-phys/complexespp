// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#ifndef ENERGY_ENERGY_H
#define ENERGY_ENERGY_H

namespace energy {

double harmonic(const double r, const double rmin, const double k) noexcept;

}  // namespace energy

#endif  // ENERGY_ENERGY_H
