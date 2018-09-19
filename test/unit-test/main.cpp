// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#include <fmt/format.h>
#include <gtest/gtest.h>
#include <omp.h>

#include "util/log.h"

int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);
  util::GlobalLog::setNumberOfThreads(omp_get_max_threads());
  util::Logger unitTestLog("foo.log");
  util::GlobalLog::setGlobalLog(&unitTestLog);
  return RUN_ALL_TESTS();
}
