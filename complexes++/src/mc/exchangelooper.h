// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
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

    bool operator!=(const Iterator& other) { return m_idx1 != other.m_idx1; }

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
}

#endif
