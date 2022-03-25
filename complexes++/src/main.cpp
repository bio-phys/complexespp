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
#include "setup/copyright.h"
#include "setup/application.h"
#include "util/log.h"
#include "util/timer.h"

int main(int argc, char **argv) {
  setup::printGPL();

  auto args = setup::CLIArgs(argc, argv);
  if (args.parse() == false || args.hasKey("help")) {
    args.printHelp(std::cout);
    return 0;
  }
  setup::Application app(args);
  auto timer = util::Timer();
  auto res = app.run();
  fmt::print(std::clog, "[LOG] total runtime = {}s\n",
             timer.stopAndGetElapsed());
  return res;
}
