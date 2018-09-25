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
#include "setup/config.h"
#include "util/file.h"

namespace setup {

Config::Config(const std::string &configFile) {
  util::throwIfFileDoesNotExists(configFile);
  m_config = YAML::LoadFile(configFile);
}

Config::Config(io::Deserializer &deserializer)
    : m_configFile(
          deserializer.restore<decltype(m_configFile)>("m_configFile")),
      m_config(m_configFile) {}

bool Config::hasValue(const std::string &key) const {
  auto const tokens = util::splitStr(key, ".");
  auto node = YAML::Clone(m_config);
  for (auto const &t : tokens) {
    node = node[t];
    if (!node.IsDefined()) {
      return false;
    }
  }
  return true;
}

void Config::serialize(io::Serializer &serializer) const {
  serializer.append(m_configFile, "m_configFile");
}
} // namespace setup
