// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#include "io/reader.h"

BaseReader::BaseReader(std::shared_ptr<domains::Domains> dom,
                       const util::rvec& box)
    : m_dom(dom), m_box(box), m_numAtomsModel(0) {
  for (auto& d : *m_dom) {
    m_numAtomsModel += d->nBeads();
  }
}

BaseReader::~BaseReader() {}
