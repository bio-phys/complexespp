// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#ifndef CONSTANTS_H
#define CONSTANTS_H

namespace constants {

namespace natural {
auto constexpr elementaryCharge = 1.602176565e-19;  // couloumb
auto constexpr epsilon_0 = 8.854187817620389e-12;   // m^-3 kg^-1 s^4 A^2
auto constexpr k = 1.3806488e-23;                   // m^2 kg s^-2 K^-1
auto constexpr refT = 300.0;                        // Kelvin
}  // namespace natural

namespace units {
auto constexpr nano = 1e-9;
auto constexpr angstrom = 1e-10;
auto constexpr energy = natural::k * natural::refT;
auto constexpr centi = 1e-2;
auto constexpr bar = 1e5;  // bar to Pa
}  // namespace units

namespace conversions {
auto constexpr kcalPerMolToJoulePerMolecule = 6.9477e-21;
auto constexpr barTOktperCubAA =
    units::bar / units::energy *
    (units::angstrom * units::angstrom * units::angstrom);
}  // namespace conversions

namespace polymerchain {
auto constexpr bondLength = 3.81;  // Angstrom
auto constexpr kPseudoBond = 378 * conversions::kcalPerMolToJoulePerMolecule /
                             units::energy;  // in kcal / mol
}  // namespace polymerchain
}  // namespace constants
#endif  // CONSTANTS_H
