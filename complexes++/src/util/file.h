// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#ifndef UTIL_FILE_H
#define UTIL_FILE_H

#include <string>
class string;

namespace util {

std::string backUpIfExists(const std::string& file);
void deleteFileIfExists(const std::string& file);
std::string basename(const std::string& file);
std::string fileSuffix(const std::string& file);
std::string absolutePath(const std::string& path);
std::string appendPaths(const std::string& basePath,
                        const std::string& subPath);
template <class... Params>
std::string appendPaths(const std::string& basePath, const std::string& subPath,
                        const std::string& morePath, Params... params) {
  return appendPaths(appendPaths(basePath, subPath), morePath, params...);
}
std::string appendPathsIfSubRelative(const std::string& basePath,
                                     const std::string& subPath);
std::string dirname(const std::string& pathWithFilename);
std::string filename(const std::string& pathWithFilename);
void throwIfFileDoesNotExists(const std::string& file);
bool isDirectory(const std::string& dirName);
}  // namespace util

#endif  // UTIL_FILE_H
