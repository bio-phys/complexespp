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
#ifndef UTIL_LOG_H
#define UTIL_LOG_H

#define FMT_HEADER_ONLY
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

#include "util/util.h"

namespace setup {
class Config;
} // namespace setup

namespace util {

class Logger {
public:
  explicit Logger() : m_stream(nullptr), m_prefix("[LOG] "), m_toCLOG(false) {}
  explicit Logger()
      : m_stream(std::make_unique<std::ofstream>(inFname, std::ofstream::app)),
        m_prefix("[LOG] "), m_toCLOG(true) {}
  explicit Logger(const std::string &inPrefix, const bool inToCLOG)
      : m_stream(),
        m_prefix(inPrefix), m_toCLOG(inToCLOG) {}

  // support moving
  Logger(Logger &&rhs) = default;
  Logger &operator=(Logger &&rhs) = default;
  // disable copying
  Logger(const Logger &rhs) = delete;
  Logger &operator=(const Logger &rhs) = delete;

  template <class... Params>
  void operator()(const std::string &inFormat, const Params &... params) {
    if (m_stream) {
      DEBUG_ASSERT(m_stream != nullptr, "Logger not initialized");
// TODO use a mutex here
#pragma omp critical(FILEOUTPUT)
      {
        fmt::print(*m_stream, m_prefix);
        fmt::print(*m_stream, inFormat, params...);
      }
    }
    else{
#pragma omp critical(FILEOUTPUT)
        buffer << m_prefix << params...;
    }
    if (m_toCLOG) {
// Standard output should be thread safe but lets consider it is not
#pragma omp critical(STDLOGOUTPUT)
      {
        fmt::print(std::clog, m_prefix);
        fmt::print(std::clog, inFormat, params...);
      }
    }
  }

  void flush() {
// TODO use a mutex here
#pragma omp critical(FILEOUTPUT)
    m_stream->flush();
  }

  bool isSet() const { return m_stream != nullptr; }

  void setLogFile(const std::string &inFname){
      DEBUG_ASSERT(!isSet(), "Logger already set");
      m_stream = std::make_unique<std::ofstream>(inFname, std::ofstream::app);
      (*m_stream) << buffer.str();
      buffer.clear();
  }

private:
  std::unique_ptr<std::ostream> m_stream = nullptr;
  std::string m_prefix;
  bool m_toCLOG;
  std::ostringstream buffer;
};

namespace GlobalLog {
void setNumberOfThreads(const int inNbThreads);
void setGlobalLog(Logger *inSimuLog);
void redirectLog(const int inThreadId);
Logger *getGlobalLog();
} // namespace GlobalLog
// Print using the current global log
template <class... Params> void Log(Params... params) {
  (*GlobalLog::getGlobalLog())(params...);
}

std::string buildInformation();

std::string configInformation(const setup::Config &config,
                              const std::string &configFile);
std::string startingString();
std::string endingString();
} // namespace util

#endif // UTIL_LOG_H
