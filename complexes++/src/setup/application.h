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
  Application(const int inArgc, char **inArgv) : m_args(inArgc, inArgv) {}

  //! Forbid copy
  Application(const Application &) = delete;
  Application &operator=(const Application &) = delete;

  //! execute the application
  int run();
};
} // namespace setup

#endif // SETUP_APPLICATION_H
