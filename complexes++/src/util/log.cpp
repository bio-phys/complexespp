// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------

#include "version.h"

#include "complexesconfig.h"
#include "parallelization/ompmanager.h"
#include "setup/config.h"
#include "util/log.h"

namespace util {

namespace GlobalLog {
std::vector<Logger*> CurrentGlobalLog(1, nullptr);

void setNumberOfThreads(const int inNbThreads) {
  CurrentGlobalLog.resize(inNbThreads, nullptr);
}

void setGlobalLog(Logger* inSimuLog) {
  CurrentGlobalLog[omp_get_thread_num()] = inSimuLog;
}

void redirectLog(const int inThreadId) {
  CurrentGlobalLog[omp_get_thread_num()] = CurrentGlobalLog[inThreadId];
}

Logger* getGlobalLog() {
  if (CurrentGlobalLog[omp_get_thread_num()] == nullptr) {
    throw std::runtime_error(
        "The global log is used but have not been initialized.");
  }
  return CurrentGlobalLog[omp_get_thread_num()];
}
}  // namespace GlobalLog

std::string buildInformation() {
  return fmt::format(
      "\n"
      "[LOG] #######################\n"
      "[LOG] # PROGRAM INFORMATION #\n"
      "[LOG] #######################\n"
      "[LOG] name: Complexes++\n"
      "[LOG] version: {}\n"
      "[LOG] commit: {}\n"
      "[LOG] branch: {}\n"
      "[LOG] compilation flags: {}\n"
      "[LOG] compilation libs: {}\n"
      "[LOG] \n",
      VERSION, GIT_COMMIT_HASH, GIT_BRANCH, ComplexesCompileFlags(),
      ComplexesCompileLibs());
}

std::string configInformation(const setup::Config& config,
                              const std::string& configFile) {
  return fmt::format(
      "\n"
      "[LOG] #################################\n"
      "[LOG] # Parsing Structure and Options #\n"
      "[LOG] #################################\n"
      "[LOG] structure file: {}\n"
      "[LOG] trajectory file: {}\n"
      "[LOG] config file: {}\n",
      config.value<std::string>("structure"),
      config.value<std::string>("output.file"), configFile);
}

std::string startingString() {
  return "\n"
         "#######################\n"
         "# Starting Simulation #\n"
         "#######################\n";
}

std::string endingString() {
  return std::string(
      "\n"
      "#####################\n"
      "# Ending Simulation #\n"
      "#####################\n");
}
}  // namespace util
