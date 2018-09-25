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
#ifndef IO_XTC_H
#define IO_XTC_H

#include "domains/abstractdomain.h"
#include "io/reader.h"

class XDRFILE;

namespace io {

// minimal wrapper around the XDRFILE library for RAII. You can give this to
// xdrfile library functions using the & operator. Smart Pointers would also be
// an option. But the code calling the xtc routines just looks nicer with this
// class.
class XDR {
public:
  explicit XDR(const std::string &file, const std::string &mode);
  ~XDR();
  // class movable
  XDR(XDR &&rhs) = default;
  XDR &operator=(XDR &&rhs) = default;
  // class not copyable
  XDR(const XDR &rhs) = delete;
  XDR &operator=(const XDR &rhs) = delete;

  XDRFILE *operator&();
  void flush();

private:
  const std::string m_file;
  XDRFILE *m_pXdrfile;
};

void writeXTC(XDR &xdr, const domains::Domains &model, const util::rvec &box,
              const int step, const double time);
void writeTRR(XDR &xdr, const domains::Domains &model, const util::rvec &box,
              const int step, const double time);

enum class XDR_TYPE { xtc, trr };

class XDRReader : public BaseReader {
public:
  explicit XDRReader(std::shared_ptr<domains::Domains> dom,
                     const util::rvec &box, const std::string &file);
  ~XDRReader();

  std::shared_ptr<domains::Domains> nextFrame() final;
  bool hasNextFrame() const final;

private:
  void readXDRFrame();

  mutable XDR m_xdr;
  XDR_TYPE m_type;
  bool m_bReadLast;
  util::rArray m_lastFrame;
};

} // namespace io

#endif // IO_XTC_H
