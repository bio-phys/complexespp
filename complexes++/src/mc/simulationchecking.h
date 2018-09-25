// Copyright (c) 2018 the complexes++ development team and contributors
// (see the file AUTHORS for the full list of names)
//
// This file is part of complexes++.
//
// complexes++ is free software: you can redistribute it and/or modify
// it under the terms of the Lesser GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// complexes++ is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with complexes++.  If not, see <https://www.gnu.org/licenses/>
#ifndef MC_SIMULATIONCHECKING_H
#define MC_SIMULATIONCHECKING_H

#include "mc/simulation.h"

namespace mc {
namespace SimulationChecking {

// easy access to enum type
template <typename E> constexpr auto toUType(E enumerator) noexcept {
  return static_cast<std::underlying_type_t<E>>(enumerator);
}

enum class ComparisonFlag {
  NONE = (1 << 0),
  NB_DOMS = (1 << 1),
  FORCEFIELD = (1 << 2),
  NB_SWEEP = (1 << 3),
};

inline ComparisonFlag operator|(ComparisonFlag lhs, ComparisonFlag rhs) {
  return ComparisonFlag(toUType(lhs) | toUType(rhs));
}

inline ComparisonFlag operator&(ComparisonFlag lhs, ComparisonFlag rhs) {
  return ComparisonFlag(toUType(lhs) & toUType(rhs));
}

inline bool operator!=(ComparisonFlag lhs, int rhs) {
  return (toUType(lhs) != rhs);
}

//! This function should be used to ensure that the given configurations
//! are compatible with the multidir exchange mode
inline bool compareSimulations(const std::vector<Simulation> &simus,
                               const ComparisonFlag comparisonMask,
                               const bool assertOnFailure) {
  if (simus.size() == 0) {
    return true;
  }

  const Simulation &firstSimu = simus[0];

  // Compare first with all others
  for (const Simulation &currentSimu : simus) {
    if ((comparisonMask & ComparisonFlag::NB_DOMS) != 0) {
      const size_t nbDomains0 = firstSimu.getNbDomains();
      const size_t nbDomainsCurrent = currentSimu.getNbDomains();
      if (nbDomains0 != nbDomainsCurrent) {
        if (assertOnFailure) {
          throw std::runtime_error(
              fmt::format("Error in multidir, simulations are not compatible.\n"
                          "First config file {} is with {} domains whereas,"
                          " config {} is with {} domains.",
                          firstSimu.getShortName(), nbDomains0,
                          currentSimu.getShortName(), nbDomainsCurrent));
        }
        return false;
      }
    }
    if ((comparisonMask & ComparisonFlag::FORCEFIELD) != 0) {
      const auto &forceField0 = firstSimu.getForceField();
      const auto &forceFieldCurrent = currentSimu.getForceField();
      if (forceField0 != forceFieldCurrent) {
        if (assertOnFailure) {
          throw std::runtime_error(fmt::format(
              "Error in multidir, simulations are not compatible.\n"
              "Force field are different between {} and {}.",
              firstSimu.getShortName(), currentSimu.getShortName()));
        }
        return false;
      }
    }
    if ((comparisonMask & ComparisonFlag::NB_SWEEP) != 0) {
      const int nbSweeps0 = firstSimu.getNbSweep();
      const int nbSweepsCurrent = currentSimu.getNbSweep();
      if (nbSweeps0 != nbSweepsCurrent) {
        if (assertOnFailure) {
          throw std::runtime_error(
              fmt::format("Error in multidir, simulations are not compatible.\n"
                          "Nb sweeps are different, for {} it is {} sweeps and "
                          "for {} it is {} sweeps.",
                          firstSimu.getShortName(), nbSweeps0,
                          currentSimu.getShortName(), nbSweepsCurrent));
        }
        return false;
      }
    }
  }
  return true;
}

inline ComparisonFlag MultiDirMask() { return ComparisonFlag::NONE; }

inline ComparisonFlag ExchangeMask() {
  return ComparisonFlag::NB_SWEEP | ComparisonFlag::NB_DOMS;
}
} // namespace SimulationChecking
} // namespace mc

#endif // MC_SIMULATIONCHECKING_H
