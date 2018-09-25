// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#include "setup/application.h"
#include "util/log.h"
#include "util/timer.h"

int main(int argc, char** argv) {
  setup::Application app(argc, argv);
  auto timer = util::Timer();
  auto res = app.run();
  fmt::print(std::clog, "[LOG] total runtime = {}s\n",
             timer.stopAndGetElapsed());
  return res;
}
