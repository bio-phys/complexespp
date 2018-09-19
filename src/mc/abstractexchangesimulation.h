// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#ifndef ABSTRACTEXCHANGEBUILDER_H
#define ABSTRACTEXCHANGEBUILDER_H

#include <string>

namespace mc {

//! This abstract class represents the interface
//! to a high level exchange alorithm.
//! The current usage is to call:
//! - init : to load from disk and restart if needed
//! - compareSimulations : to ensure compatible configs
//! - run : to execute the complete exchange algo
class AbstractExchangeSimulation {
 public:
  virtual ~AbstractExchangeSimulation() {}
  virtual bool compareSimulations(const bool assertOnFailure) = 0;
  virtual int run() = 0;
};
}

#endif
