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
#ifndef IO_STATFILE_H
#define IO_STATFILE_H

#include <fstream>
#include <string>

#include "io/serializer.h"

namespace io {

class StatFile : public AbstractSerializable {
public:
  StatFile(const std::string fname)
      : m_fname(fname), m_out(m_fname.c_str(), std::ios_base::out) {}

  // class movable
  StatFile(StatFile &&rhs) = default;
  StatFile &operator=(StatFile &&rhs) = default;
  // class not copyable
  StatFile(const StatFile &rhs) = delete;
  StatFile &operator=(const StatFile &rhs) = delete;

  std::ostream &stream() { return m_out; }
  void flush() { m_out.flush(); }
  const std::string &fname() const noexcept { return m_fname; };

  void serialize(io::Serializer &serializer) const final {
    serializer.append(m_fname, "m_fname");
  };

  StatFile(io::Deserializer &deserializer)
      : m_fname(deserializer.restore<decltype(m_fname)>("m_fname")),
        m_out(m_fname, std::ios_base::app) {}

private:
  std::string m_fname;
  std::ofstream m_out;
};

// Base case
template <class T>
typename std::enable_if<util::is_std_vector<T>::value, void>::type
writeStats(StatFile &out, const T &arg) {
  for (unsigned i = 0; i + 1 < arg.size(); ++i) {
    fmt::print(out.stream(), "{},", arg[i]);
  }
  if (arg.size()) {
    fmt::print(out.stream(), "{}", arg[arg.size() - 1]);
  }
  fmt::print(out.stream(), "\n");
}

template <class T>
typename std::enable_if<!util::is_std_vector<T>::value, void>::type
writeStats(StatFile &out, const T &arg) {
  fmt::print(out.stream(), "{}\n", arg);
}

template <class T, class... Params>
typename std::enable_if<util::is_std_vector<T>::value, void>::type
writeStats(StatFile &out, const T &head, const Params &... tail) {
  for (unsigned i = 0; i < head.size(); ++i) {
    fmt::print(out.stream(), "{},", head[i]);
  }
  writeStats(out, tail...);
}

template <class T, class... Params>
typename std::enable_if<!util::is_std_vector<T>::value, void>::type
writeStats(StatFile &out, const T &head, const Params &... tail) {
  fmt::print(out.stream(), "{},", head);
  writeStats(out, tail...);
}

template <class... Params>
void writeStatsHeader(StatFile &out, const Params &... params) {
  fmt::print(out.stream(), "# ");
  writeStats(out, params...);
}

} // namespace io

#endif // IO_STATFILE_H
