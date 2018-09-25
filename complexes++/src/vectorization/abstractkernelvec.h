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
#ifndef ABSTRACTKERNELVEC_H
#define ABSTRACTKERNELVEC_H

namespace simd {

template <class VecType> class AbstractKernelVec {
public:
  /** Return the energy using Lennard Jones and Hueckel formulation
   * Does not have to be const or thread safe.
   */
  virtual VecType compute(const VecType r2, const VecType inCharge1,
                          const VecType inCharge2, const int inBeadType1,
                          const int inBeadType2[]) = 0;
};
} // namespace simd

#endif
