// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#ifndef REBUILDMEMBRANE_HPP
#define REBUILDMEMBRANE_HPP

#include "membrane/abstractmembrane.h"
#include "membrane/flatmembrane.h"
#include "membrane/spheremembrane.h"
#include "membrane/tubemembrane.h"

namespace membrane {

template <typename T>
static std::unique_ptr<membrane::AbstractMembrane<T>> RebuildMembrane(
    io::Deserializer& deserializer, const std::string key) {
  deserializer.access(key);
  return AbstractMembrane<T>::Rebuild(deserializer);
}
}  // namespace membrane

#endif
