// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#ifndef EXCHANGEACCEPTER_H
#define EXCHANGEACCEPTER_H

#include "mc/simulation.h"

namespace mc {

//! This file contains different accept classes to be used
//! by the exchange algorithm.
//! Each class must propose a method "accept" that takes in parameters:
//! - const int inIdx1, the index of the first simulation
//! - const int inIdx2, the index of the second simulation
//! - const setup::Simulation& inS1, the first simulation class
//! - const setup::Simulation& inS2, the second simulation class
//! - util::RNGEngine& inOutRng, a random generator to be used if the class
//! needs
//! some
//! - double* outProbability, a pointer to store the probability
//!
//! In return the method "accept" can return one of the two options:
//! - std::tuple<bool>, where the boolean tell if the exchange is accepted
//! - std::tuple<bool,util::rArray,util::rArray>, where the boolean
//! tell if the exchange is accepted, and the two arrays are the energy for each
//! simulation using the forcefield of the other
//! This second option should be used when the "accept" method have computed the
//! two arrays
//! and wants to return the result to potentially avoid recomputation
//! If the accept if "false" the energy matrices returned could be empty
//! (it is accept function implementation choice).

/**
 * See:
 * Bussi, G. (2014). Hamiltonian replica exchange in GROMACS: a flexible
 * implementation.
 * Molecular Physics, 112(3–4), 379–384.
 * http://doi.org/10.1080/00268976.2013.824126
 * for more details.
 */
class HREXAccept {
 public:
  bool needSameForceField() const { return false; }

  std::tuple<bool, double, energy::rEnergyMatrix, energy::rEnergyMatrix> accept(
      const Simulation& inS1, const Simulation& inS2,
      util::RNGEngine& inOutRng) {
    const double Ui_ri = inS1.getEnergy();
    energy::rEnergyMatrix Ui_rj_array =
        inS2.computeEnergyForFF(inS1.getForceField());
    const double Ui_rj = Ui_rj_array.getTotalEnergy();
    const double Uj_rj = inS2.getEnergy();
    energy::rEnergyMatrix Uj_ri_array =
        inS1.computeEnergyForFF(inS2.getForceField());
    const double Uj_ri = Uj_ri_array.getTotalEnergy();
    const double Bi = inS1.getBeta();
    const double Bj = inS2.getBeta();
    const double coef =
        std::min(1., exp(((Ui_ri - Ui_rj) * Bi) + ((Uj_rj - Uj_ri) * Bj)));
    return std::make_tuple(
        std::uniform_real_distribution<>{0, 1}(inOutRng) < coef, coef,
        std::move(Uj_ri_array), std::move(Ui_rj_array));
  }
};

/**
 * See :
 * Hansmann, U. H. E. (1997).
 * Parallel tempering algorithm for conformational studies of biological
 * molecules.
 * Chem. Phys. Lett., 281(1–3), 140–150.
 * http://doi.org/10.1016/S0009-2614(97)01198-6
 * for more details.
 */
class REMCAccept {
 public:
  bool needSameForceField() const { return true; }

  std::tuple<bool, double> accept(const Simulation& inS1,
                                  const Simulation& inS2,
                                  util::RNGEngine& inOutRng) {
    const double deltaE = inS2.getEnergy() - inS1.getEnergy();
    const double B = inS2.getBeta() - inS1.getBeta();
    const double coef = std::min(1., exp(B * deltaE));
    return std::make_tuple(
        std::uniform_real_distribution<>{0, 1}(inOutRng) < coef, coef);
  }
};

/**
 * See :
 * Okabe, T., Kawata, M., Okamoto, Y., & Mikami, M. (2001). Replica-exchange
 * Monte Carlo method for the isobaric-isothermal ensemble. Chemical Physics
 * Letters, 335(5–6), 435–439. https://doi.org/10.1016/S0009-2614(01)00055-0
 * for more details.
 */
class NPTAccept {
 public:
  bool needSameForceField() const { return true; }

  std::tuple<bool, double> accept(const Simulation& inS1,
                                  const Simulation& inS2,
                                  util::RNGEngine& inOutRng) {
    const double deltaE = inS1.getEnergy() - inS2.getEnergy();
    const double B = inS1.getBeta() - inS2.getBeta();
    const auto dV = util::volume(inS1.getBox()) - util::volume(inS2.getBox());
    const auto dP = inS1.getBeta() * inS1.getPressure() -
                    inS2.getBeta() * inS2.getPressure();
    const double coef = std::min(1., exp(B * deltaE + dP * dV));
    return std::make_tuple(
        std::uniform_real_distribution<>{0, 1}(inOutRng) < coef, coef);
  }
};

class TrueAccept {
 public:
  bool needSameForceField() const { return false; }

  std::tuple<bool, double> accept(const Simulation& /*inS1*/,
                                  const Simulation& /*inS2*/,
                                  util::RNGEngine& /*inOutRng*/) {
    return std::make_tuple(true, 1);
  }
};

class FalseAccept {
 public:
  bool needSameForceField() const { return false; }

  std::tuple<bool, double> accept(const Simulation& /*inS1*/,
                                  const Simulation& /*inS2*/,
                                  util::RNGEngine& /*inOutRng*/) {
    return std::make_tuple(false, 0);
  }
};
}  // namespace mc

#endif
