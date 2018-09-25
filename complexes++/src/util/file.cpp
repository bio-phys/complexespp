// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#include <boost/filesystem.hpp>
#include <fmt/format.h>
#include <fmt/time.h>
#include <string>

#include "util/file.h"

namespace fs = boost::filesystem;

namespace util {

// rename file to backup-POSIXTIME-originalname if the file already exists.
std::string backUpIfExists(const std::string& file) {
  if (!fs::exists(file)) {
    return std::string();
  }
  fs::path bpath(file);
  auto datetime = std::time(nullptr);
  auto newName = appendPaths(
      bpath.parent_path().string(),
      fmt::format("backup_{0:%Y-%m-%d}:{0:%H-%M-%S}_{1}",
                  *std::localtime(&datetime), bpath.filename().string()));
  fs::rename(file, newName);
  return newName;
}

void deleteFileIfExists(const std::string& file) {
  if (fs::exists(file)) {
    fs::remove(file);
  }
}

std::string basename(const std::string& file) {
  return fs::basename(file);
}

void throwIfFileDoesNotExists(const std::string& file) {
  if (!fs::exists(file)) {
    throw std::invalid_argument(fmt::format("{} <-- does not exist.\n", file));
  }
}

std::string fileSuffix(const std::string& file) {
  return fs::path(file).extension().string();
}

std::string absolutePath(const std::string& path) {
  return fs::canonical(fs::path(path), fs::current_path()).string();
}

std::string appendPaths(const std::string& basePath,
                        const std::string& subPath) {
  return (fs::path(basePath) / fs::path(subPath)).string();
}

std::string appendPathsIfSubRelative(const std::string& basePath,
                                     const std::string& subPath) {
  fs::path bfsSubPath(subPath);
  if (bfsSubPath.is_relative()) {
    return (fs::path(basePath) / bfsSubPath).string();
  } else {
    return subPath;
  }
}

std::string filename(const std::string& pathWithFilename) {
  return fs::path(pathWithFilename).filename().string();
}

std::string dirname(const std::string& pathWithFilename) {
  fs::path file(pathWithFilename);
  if (file.has_filename()) {
    return file.parent_path().string();
  }
  return pathWithFilename;
}

bool isDirectory(const std::string& dirName) {
  return fs::is_directory(dirName);
}
}  // namespace
