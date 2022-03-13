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
#ifndef SETUP_CLIARGS_H
#define SETUP_CLIARGS_H

#define FMT_HEADER_ONLY
#include <any>
#include <fmt/format.h>
#include <initializer_list>
#include <utility>
#include <unordered_map>

class string;
class ostream;

namespace setup {
class CLIArgs {
protected:
  CLIArgs();

public:
  CLIArgs(const int &argc, const char *const argv[]);

  template <class T> T value(const std::string &key) const {
    auto val = T();
    try {
      val = std::any_cast<T>(m_args.at(key));
    } catch (...) {
      throw std::invalid_argument(
          fmt::format("No CLI argument called : {}", key));
    }

    return val;
  }
  std::string value(const std::string &key) const;
  bool hasKey(const std::string &key) const noexcept;

  template <class T>
  T getMappingValue(
      const std::string &key,
      std::initializer_list<std::pair<std::string, T>> inMapping) const {
    const std::string strValue = value(key);

    for (const auto &keyvalue : inMapping) {
      if (keyvalue.first == strValue) {
        return keyvalue.second;
      }
    }

    // Build list for error
    std::string allValidValues;
    for (const auto &keyvalue : inMapping) {
      allValidValues += keyvalue.first + ", ";
    }

    throw std::invalid_argument(
        fmt::format("Invalid argument for parameter {}, you passed {} but only "
                    "{} are valid ",
                    key, strValue, allValidValues));
  }

protected:
  std::unordered_map<std::string, std::any> m_args;
};

} // namespace setup

#endif // SETUP_CLIARGS_H
