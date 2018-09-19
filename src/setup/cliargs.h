// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#ifndef SETUP_CLIARGS_H
#define SETUP_CLIARGS_H

#include <boost/program_options.hpp>
#include <fmt/format.h>
#include <initializer_list>
#include <utility>

class string;
class ostream;

namespace setup {
class CLIArgs {
 protected:
  CLIArgs();

 public:
  CLIArgs(const int& argc, const char* const argv[]);

  template <class T>
  T value(const std::string& key) const {
    auto val = T();
    try {
      val = m_args[key].as<T>();
    } catch (...) {
      throw std::invalid_argument(
          fmt::format("No CLI argument called : {}", key));
    }

    return val;
  }
  std::string value(const std::string& key) const;
  bool hasKey(const std::string& key) const noexcept;
  std::ostream& print(std::ostream& out) const;

  template <class T>
  T getMappingValue(
      const std::string& key,
      std::initializer_list<std::pair<std::string, T>> inMapping) const {
    const std::string strValue = value(key);

    for (const auto& keyvalue : inMapping) {
      if (keyvalue.first == strValue) {
        return keyvalue.second;
      }
    }

    // Build list for error
    std::string allValidValues;
    for (const auto& keyvalue : inMapping) {
      allValidValues += keyvalue.first + ", ";
    }

    throw std::invalid_argument(
        fmt::format("Invalid argument for parameter {}, you passed {} but only "
                    "{} are valid ",
                    key, strValue, allValidValues));
  }

 protected:
  void defineArgs();

  boost::program_options::options_description m_required;
  boost::program_options::variables_map m_args;
};

void printHelp(const CLIArgs& args);

std::ostream& operator<<(std::ostream& s, const CLIArgs&);
}  // namespace setup

#endif  // SETUP_CLIARGS_H
