// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#include <boost/algorithm/string.hpp>
#include <string>
#include <vector>

#include "util/log.h"
#include "util/string.h"

namespace util {

std::vector<std::string> splitStr(const std::string& str,
                                  const std::string& delimiter) {
  std::vector<std::string> tokens;
  const auto trimmed = boost::trim_copy(str);
  boost::split(tokens, trimmed, boost::is_any_of(delimiter),
               boost::algorithm::token_compress_on);
  return tokens;
}

int truncateLeft(const int val, const int ndigits) {
  const int base = 10;
  const int max = static_cast<int>(std::pow(base, ndigits)) - 1;
  if (val <= max) {
    return val;
  }
  int mul = 1, leftDigit = val;
  while (leftDigit >= base) {
    leftDigit /= base;
    mul *= base;
  }
  return val - (leftDigit * mul);
}
}  // namespace util
