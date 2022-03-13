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
#include <filesystem>
#define FMT_HEADER_ONLY
#include <fmt/format.h>
#include <fmt/chrono.h>
#include <string>

#include "util/file.h"

namespace fs = std::filesystem;

namespace util {

// rename file to backup-POSIXTIME-originalname if the file already exists.
std::string backUpIfExists(const std::string &file) {
  if (!fs::exists(file)) {
    return std::string();
  }
  fs::path bpath(file);
  auto datetime = std::time(nullptr);
  auto newName = appendPaths(bpath.parent_path().string(),
                             fmt::format("backup_{0:%Y-%m-%d}:{0:%H-%M-%S}_{1}",
                                         *std::localtime(&datetime),
                                         bpath.filename().string()));
  fs::rename(file, newName);
  return newName;
}

void deleteFileIfExists(const std::string &file) {
  if (fs::exists(file)) {
    fs::remove(file);
  }
}

std::string basename(const std::string &file) { return fs::path(file).stem().string(); }

void throwIfFileDoesNotExists(const std::string &file) {
  if (!fs::exists(file)) {
    throw std::invalid_argument(fmt::format("{} <-- does not exist.\n", file));
  }
}

std::string fileSuffix(const std::string &file) {
  return fs::path(file).extension().string();
}

std::string absolutePath(const std::string &path) {
  return fs::absolute(fs::path(path)).string();
}

std::string appendPaths(const std::string &basePath,
                        const std::string &subPath) {
  return (fs::path(basePath) / fs::path(subPath)).string();
}

std::string appendPathsIfSubRelative(const std::string &basePath,
                                     const std::string &subPath) {
  fs::path bfsSubPath(subPath);
  if (bfsSubPath.is_relative()) {
    return (fs::path(basePath) / bfsSubPath).string();
  } else {
    return subPath;
  }
}

std::string filename(const std::string &pathWithFilename) {
  return fs::path(pathWithFilename).filename().string();
}

std::string dirname(const std::string &pathWithFilename) {
  fs::path file(pathWithFilename);
  if (file.has_filename()) {
    return file.parent_path().string();
  }
  return pathWithFilename;
}

bool isDirectory(const std::string &dirName) {
  return fs::is_directory(dirName);
}
} // namespace util
