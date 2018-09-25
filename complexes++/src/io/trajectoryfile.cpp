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
#include "io/trajectoryfile.h"
#include "io/serializer.h"
#include "util/file.h"

namespace io {
std::string deduceType(const std::string &fname) {
  auto suffix = util::fileSuffix(fname);
  suffix.erase(0, 1);
  return suffix;
}

std::string modeToXDR(const TrajectoryFile::Mode mode) {
  switch (mode) {
  case TrajectoryFile::Mode::app:
    return "a";
  case TrajectoryFile::Mode::in:
    return "r";
  case TrajectoryFile::Mode::out:
    return "w";
  }
  throw std::runtime_error("unreachable code to quite gcc warnings");
}

std::ios_base::openmode modeToFStream(const TrajectoryFile::Mode mode) {
  switch (mode) {
  case TrajectoryFile::Mode::app:
    return std::ios_base::app;
  case TrajectoryFile::Mode::in:
    return std::ios_base::in;
  case TrajectoryFile::Mode::out:
    return std::ios_base::out;
  }
  throw std::runtime_error("unreachable code to quite gcc warnings");
}

TrajectoryFile::TrajectoryFile(const std::string &_fname,
                               const TrajectoryFile::Mode _mode)
    : m_fname(_fname), m_type(deduceType(_fname)), m_mode(_mode) {
  if (m_type == "xtc" || m_type == "trr") {
    m_xdr = std::make_unique<XDR>(m_fname, modeToXDR(m_mode));
  } else if (m_type == "pdb" || m_type == "gro") {
    m_fstream = std::make_unique<std::ofstream>(m_fname, modeToFStream(m_mode));
  } else {
    throw std::invalid_argument(
        fmt::format("don't know that trajectory format: {}", m_type));
  }
}

void TrajectoryFile::serialize(io::Serializer &serializer) const {
  serializer.append(m_fname, "m_fname");
  serializer.append(m_type, "m_type");
  serializer.append(m_mode, "m_mode");
}

TrajectoryFile::TrajectoryFile(io::Deserializer &deserializer)
    : m_fname(deserializer.restore<decltype(m_fname)>("m_fname")),
      m_type(deserializer.restore<decltype(m_type)>("m_type")),
      m_mode(deserializer.restore<decltype(m_mode)>("m_mode")) {
  // append to files that are written. Other modes do not need to be changed
  // on a restore
  if (m_mode == TrajectoryFile::Mode::out) {
    m_mode = TrajectoryFile::Mode::app;
  }
  if (m_mode == TrajectoryFile::Mode::in) {
    throw std::runtime_error(
        "restart with reading a trajectory file is not supported");
  }
  if (m_type == "xtc" || m_type == "trr") {
    m_xdr = std::make_unique<XDR>(m_fname, modeToXDR(m_mode));
  } else if (m_type == "pdb" || m_type == "gro") {
    m_fstream = std::make_unique<std::ofstream>(m_fname, modeToFStream(m_mode));
  } else {
    throw std::invalid_argument(
        fmt::format("don't know that trajectory format: {}", m_type));
  }
}

XDR &TrajectoryFile::xdr() {
  if (m_xdr == nullptr) {
    throw std::invalid_argument(
        fmt::format("TrajectoryFile doesn't contain a xdr"));
  }
  return *m_xdr;
}

std::ofstream &TrajectoryFile::fstream() {
  if (m_fstream == nullptr) {
    throw std::invalid_argument(
        fmt::format("TrajectoryFile doesn't contain a stream"));
  }
  return *m_fstream;
}

const std::string &TrajectoryFile::fname() const noexcept { return m_fname; }

const std::string &TrajectoryFile::type() const noexcept { return m_type; }

void TrajectoryFile::flush() {
  if (m_fstream != nullptr) {
    m_fstream->flush();
  } else if (m_xdr != nullptr) {
    m_xdr->flush();
  }
}
} // namespace io
