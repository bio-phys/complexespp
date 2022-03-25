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
#ifndef COPYRIGHT_H
#define COPYRIGHT_H

#include <iostream>

namespace setup {

inline void printGPL() {
  std::cout << " Copyright (c) 2018 the complexes++ development team and "
               "contributors \n";
  std::cout << " (see the file AUTHORS for the full list of names) \n";
  std::cout << " \n";
  std::cout << " complexes++ is free software: you can redistribute it and/or "
               "modify \n";
  std::cout << " it under the terms of the Lesser GNU General Public License "
               "as published by \n";
  std::cout << " the Free Software Foundation, either version 3 of the "
               "License, or \n";
  std::cout << " (at your option) any later version. \n";
  std::cout << " \n";
  std::cout
      << " complexes++ is distributed in the hope that it will be useful, \n";
  std::cout
      << " but WITHOUT ANY WARRANTY; without even the implied warranty of \n";
  std::cout
      << " MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the \n";
  std::cout << " GNU General Public License for more details. \n";
  std::cout << " \n";
  std::cout << " You should have received a copy of the GNU General Public "
               "License \n";
  std::cout << " along with complexes++.  If not, see "
               "<https://www.gnu.org/licenses/> \n";
  std::cout << std::endl;
}

}

#endif // COPYRIGHT_H
