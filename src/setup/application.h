// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#ifndef SETUP_APPLICATION_H
#define SETUP_APPLICATION_H

#include <memory>

#include "setup/cliargs.h"

namespace setup {

/**
 * @brief The Application class is in charge of the high level
 * choice of the algorithm. It decides what kind of algorithm should
 * be executed.
 * @code Application app(argc, argv);
 * @code return app.run();
 */
class Application {
  //////////////////////////////////////////////////////////////////////////
  /// Attributes
  //////////////////////////////////////////////////////////////////////////

  const CLIArgs m_args;

  //////////////////////////////////////////////////////////////////////////
  /// Private methods
  //////////////////////////////////////////////////////////////////////////

  //! For the remc/exchange executions
  int multiDirExecution() const;
  //! For the single execution
  int singleSrcExecution() const;

  //////////////////////////////////////////////////////////////////////////
  /// Public methods
  //////////////////////////////////////////////////////////////////////////
 public:
  Application(const int inArgc, char** inArgv) : m_args(inArgc, inArgv) {}

  //! Forbid copy
  Application(const Application&) = delete;
  Application& operator=(const Application&) = delete;

  //! execute the application
  int run();
};
}  // namespace setup

#endif  // SETUP_APPLICATION_H
