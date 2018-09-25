// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#ifndef UTIL_STRING_H
#define UTIL_STRING_H

class string;
class vector;

namespace util {

std::vector<std::string> splitStr(const std::string& str,
                                  const std::string& delimiter);

/* truncate left most number from an integer. So 12345 -> 2345
 * Useful to prevent overflow in writing a PDB/GRO
 */
int truncateLeft(const int val, const int ndigits = 5);

}  // namespace util

#endif  // UTIL_STRING_H
