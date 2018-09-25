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
#ifndef IO_TRAJECTORYFILE_H
#define IO_TRAJECTORYFILE_H

#include <fstream>

#include "io/serializer.h"
#include "io/xdr.h"

namespace io {

class Serialize;
class Deserializer;

std::string deduceType(const std::string &fname);

class TrajectoryFile : public AbstractSerializable {
public:
  enum class Mode { app, in, out };

  explicit TrajectoryFile(const std::string &fname,
                          const Mode mode = Mode::app);

  XDR &xdr();
  std::ofstream &fstream();
  const std::string &fname() const noexcept;
  const std::string &type() const noexcept;
  void flush();

  // class movable
  TrajectoryFile(TrajectoryFile &&rhs) = default;
  TrajectoryFile &operator=(TrajectoryFile &&rhs) = default;
  // class not copyable
  TrajectoryFile(const TrajectoryFile &rhs) = delete;
  TrajectoryFile &operator=(const TrajectoryFile &rhs) = delete;

  void serialize(io::Serializer &serializer) const final;
  TrajectoryFile(io::Deserializer &deserializer);

private:
  std::string m_fname;
  std::string m_type;
  Mode m_mode;
  std::unique_ptr<std::ofstream> m_fstream = nullptr;
  std::unique_ptr<XDR> m_xdr = nullptr;
};
} // namespace io
#endif // IO_TRAJECTORYFILE_H
