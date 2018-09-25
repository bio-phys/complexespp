// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#ifndef MC_ACCEPT_FUNC_H
#define MC_ACCEPT_FUNC_H

namespace mc {
double glauberAccept(const double deltaE, const double beta) noexcept;

double metropolisAccept(const double deltaE, const double beta) noexcept;

double dynamicAccept(const double deltaE, const double beta) noexcept;

double alwaysAccept(const double deltaE, const double beta) noexcept;
}  // namespace mc

#endif  // MC_ACCEPT_FUNC_H
