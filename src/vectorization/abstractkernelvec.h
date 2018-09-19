// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#ifndef ABSTRACTKERNELVEC_H
#define ABSTRACTKERNELVEC_H

namespace simd {

template <class VecType>
class AbstractKernelVec {
 public:
  /** Return the energy using Lennard Jones and Hueckel formulation
    * Does not have to be const or thread safe.
    */
  virtual VecType compute(const VecType r2, const VecType inCharge1,
                          const VecType inCharge2, const int inBeadType1,
                          const int inBeadType2[]) = 0;
};
}

#endif
