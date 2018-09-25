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
#ifndef IO_READER_H
#define IO_READER_H

#include <memory>

#include "domains/abstractdomain.h"
#include "util/array.h"

class BaseReader {
public:
  explicit BaseReader(std::shared_ptr<domains::Domains> dom,
                      const util::rvec &box);
  virtual ~BaseReader();
  // support moves
  BaseReader(BaseReader &&rhs) = default;
  BaseReader &operator=(BaseReader &&rhs) = default;
  // don't support copying
  BaseReader(const BaseReader &rhs) = delete;
  BaseReader &operator=(const BaseReader &rhs) = delete;

  virtual std::shared_ptr<domains::Domains> nextFrame() = 0;
  virtual bool hasNextFrame() const = 0;

protected:
  std::shared_ptr<domains::Domains> m_dom;
  const util::rvec &m_box;
  int m_numAtomsModel;
};

using Reader = std::unique_ptr<BaseReader>;

#endif // IO_READER_H
