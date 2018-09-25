// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#ifndef UTIL_LOG_H
#define UTIL_LOG_H

#include <fmt/format.h>
#include <fmt/ostream.h>
#include <fstream>
#include <iostream>
#include <string>

#include "util/util.h"

namespace setup {
class Config;
}  // namespace setup

namespace util {

class Logger {
 public:
  explicit Logger() : m_stream(nullptr), m_prefix("[LOG] "), m_toCLOG(false) {}
  explicit Logger(const std::string& inFname)
      : m_stream(std::make_unique<std::ofstream>(inFname, std::ofstream::app)),
        m_prefix("[LOG] "),
        m_toCLOG(true) {}
  explicit Logger(const std::string& inFname, const std::string& inPrefix,
                  const bool inToCLOG)
      : m_stream(std::make_unique<std::ofstream>(inFname, std::ofstream::app)),
        m_prefix(inPrefix),
        m_toCLOG(inToCLOG) {}

  // support moving
  Logger(Logger&& rhs) = default;
  Logger& operator=(Logger&& rhs) = default;
  // disable copying
  Logger(const Logger& rhs) = delete;
  Logger& operator=(const Logger& rhs) = delete;

  template <class... Params>
  void operator()(const std::string& inFormat, const Params&... params) {
    if (m_stream) {
      DEBUG_ASSERT(m_stream != nullptr, "Logger not initialized");
// TODO use a mutex here
#pragma omp critical(FILEOUTPUT)
      {
        fmt::print(*m_stream, m_prefix);
        fmt::print(*m_stream, inFormat, params...);
      }
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

 private:
  std::unique_ptr<std::ostream> m_stream = nullptr;
  std::string m_prefix;
  bool m_toCLOG;
};

namespace GlobalLog {
void setNumberOfThreads(const int inNbThreads);
void setGlobalLog(Logger* inSimuLog);
void redirectLog(const int inThreadId);
Logger* getGlobalLog();
}  // namespace GlobalLog
// Print using the current global log
template <class... Params>
void Log(Params... params) {
  (*GlobalLog::getGlobalLog())(params...);
}

std::string buildInformation();

std::string configInformation(const setup::Config& config,
                              const std::string& configFile);
std::string startingString();
std::string endingString();
}  // namespace util

#endif  // UTIL_LOG_H
