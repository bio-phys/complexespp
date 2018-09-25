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
#ifndef UTIL_FILE_H
#define UTIL_FILE_H

#include <string>
class string;

namespace util {

std::string backUpIfExists(const std::string &file);
void deleteFileIfExists(const std::string &file);
std::string basename(const std::string &file);
std::string fileSuffix(const std::string &file);
std::string absolutePath(const std::string &path);
std::string appendPaths(const std::string &basePath,
                        const std::string &subPath);
template <class... Params>
std::string appendPaths(const std::string &basePath, const std::string &subPath,
                        const std::string &morePath, Params... params) {
  return appendPaths(appendPaths(basePath, subPath), morePath, params...);
}
std::string appendPathsIfSubRelative(const std::string &basePath,
                                     const std::string &subPath);
std::string dirname(const std::string &pathWithFilename);
std::string filename(const std::string &pathWithFilename);
void throwIfFileDoesNotExists(const std::string &file);
bool isDirectory(const std::string &dirName);
} // namespace util

#endif // UTIL_FILE_H
