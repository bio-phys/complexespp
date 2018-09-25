// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#include "setup/config.h"
#include "util/file.h"

namespace setup {

Config::Config(const std::string& configFile) {
  util::throwIfFileDoesNotExists(configFile);
  m_config = YAML::LoadFile(configFile);
}

Config::Config(io::Deserializer& deserializer)
    : m_configFile(
          deserializer.restore<decltype(m_configFile)>("m_configFile")),
      m_config(m_configFile) {}

bool Config::hasValue(const std::string& key) const {
  auto const tokens = util::splitStr(key, ".");
  auto node = YAML::Clone(m_config);
  for (auto const& t : tokens) {
    node = node[t];
    if (!node.IsDefined()) {
      return false;
    }
  }
  return true;
}

void Config::serialize(io::Serializer& serializer) const {
  serializer.append(m_configFile, "m_configFile");
}
}  // namespace setup
