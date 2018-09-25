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
#ifndef EXCHANGELOOP_H
#define EXCHANGELOOP_H

#include <utility>

namespace mc {

class OddEvenNeighborLoop {
  int m_nbReplicas;
  int m_idxExchange;

public:
  OddEvenNeighborLoop(const int nbReplicas, const int idxExchange)
      : m_nbReplicas(nbReplicas), m_idxExchange(idxExchange) {}

  class Iterator {
    int m_idx1;

  public:
    Iterator(const int inIdx1) : m_idx1(inIdx1) {}

    bool operator!=(const Iterator &other) { return m_idx1 != other.m_idx1; }

    Iterator operator++() {
      m_idx1 += 2;
      return *this;
    }

    const std::pair<int, int> operator*() const { return {m_idx1, m_idx1 + 1}; }
  };

  Iterator begin() const { return Iterator(m_idxExchange & 1); }

  Iterator end() const {
    if ((m_idxExchange & 1) == (m_nbReplicas & 1)) {
      return Iterator(m_nbReplicas);
    } else {
      return Iterator(m_nbReplicas - 1);
    }
  }

  static bool IsOddEvenLoop() { return true; }
};
} // namespace mc

#endif
