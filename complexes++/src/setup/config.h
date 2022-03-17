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
#ifndef SETUP_CONFIG_H
#define SETUP_CONFIG_H

#define FMT_HEADER_ONLY
#include <fmt/format.h>
#include <yaml-cpp/yaml.h>

#include "io/serializer.h"
#include "util/log.h"
#include "util/string.h"
#include "util/demangle.h"

namespace setup {

class Config : public io::AbstractSerializable {
public:
  explicit Config(const std::string &configFile);
  explicit Config(io::Deserializer &deserializer);

  void serialize(io::Serializer &serializer) const final;

  template <class T> T value(const std::string &key) const {
    auto const tokens = util::splitStr(key, '.');

    auto node = YAML::Clone(m_config);
    for (auto const &t : tokens) {
      node = node[t];
      if (!node.IsDefined()) {
        throw std::invalid_argument(
            fmt::format("key '{}' does not exists in configuration file", key));
      }
    }

    auto val = T();
    try {
      val = node.as<T>();
    } catch (YAML::TypedBadConversion<T> &e) {
      throw std::invalid_argument(
          fmt::format("key ({}) is not of type: {}", key,
                      util::demangle(typeid(T).name())));
    }

    return val;
  }

  // Here we can store configuration values in the 'experimental' sub section of
  // the config file. These are the only values where we allow defaults as they
  // need to be turned on explicitly.
  template <typename T>
  T experimental_value(const std::string &key, const T def) const {
    auto const tokens = util::splitStr(key, '.');

    auto node = YAML::Clone(m_config)["experimental"];
    for (auto const &t : tokens) {
      node = node[t];
      if (!node.IsDefined()) {
        return def;
      }
    }

    auto val = T();
    try {
      val = node.as<T>();
    } catch (YAML::TypedBadConversion<T> &e) {
      throw std::invalid_argument(
          fmt::format("Experimental key ({}) is not of type: {}", key,
                      util::demangle(typeid(T).name())));
    }
    util::Log("Using experimental key: {}={}\n", key, val);

    return val;
  }

  bool hasValue(const std::string &key) const;

private:
  const std::string m_configFile;
  YAML::Node m_config;
};

} // namespace setup

#endif // SETUP_CONFIG_H
