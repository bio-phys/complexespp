// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#ifndef IO_READER_H
#define IO_READER_H

#include <memory>

#include "domains/abstractdomain.h"
#include "util/array.h"

class BaseReader {
 public:
  explicit BaseReader(std::shared_ptr<domains::Domains> dom,
                      const util::rvec& box);
  virtual ~BaseReader();
  // support moves
  BaseReader(BaseReader&& rhs) = default;
  BaseReader& operator=(BaseReader&& rhs) = default;
  // don't support copying
  BaseReader(const BaseReader& rhs) = delete;
  BaseReader& operator=(const BaseReader& rhs) = delete;

  virtual std::shared_ptr<domains::Domains> nextFrame() = 0;
  virtual bool hasNextFrame() const = 0;

 protected:
  std::shared_ptr<domains::Domains> m_dom;
  const util::rvec& m_box;
  int m_numAtomsModel;
};

using Reader = std::unique_ptr<BaseReader>;

#endif  // IO_READER_H
