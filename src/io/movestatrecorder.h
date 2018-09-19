// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#ifndef MOVESTATRECORDER_H
#define MOVESTATRECORDER_H

#include <sstream>
#include <unordered_map>
#include <vector>

namespace io {

class MoveStatRecorder {
 protected:
  class Events {
   protected:
    int m_nbFailMoves;
    int m_nbSucceedMoves;
    int m_nbAcceptedMoves;

   public:
    Events() : m_nbFailMoves(0), m_nbSucceedMoves(0), m_nbAcceptedMoves(0) {}

    void addEvent(const bool moveSucceed, const bool moveAccepted) {
      if (moveSucceed) {
        m_nbSucceedMoves += 1;
        if (moveAccepted) {
          m_nbAcceptedMoves += 1;
        }
      } else {
        m_nbFailMoves += 1;
      }
    }

    int getNbTotalMoves() const { return m_nbFailMoves + m_nbSucceedMoves; }

    int getNbFailMoves() const { return m_nbFailMoves; }

    int getNbSucceedMoves() const { return m_nbSucceedMoves; }

    int getNbAcceptedMoves() const { return m_nbAcceptedMoves; }

    double getNbFailMovesPerc() const {
      return 100. * double(m_nbFailMoves) / double(getNbTotalMoves());
    }

    int getNbSucceedMovesPerc() const {
      return 100. * double(m_nbSucceedMoves) / double(getNbTotalMoves());
      return m_nbSucceedMoves;
    }

    int getNbAcceptedMovesPerc() const {
      return 100. * double(m_nbAcceptedMoves) / double(getNbTotalMoves());
    }
  };

  std::unordered_map<std::string, Events> m_movePerType;
  std::vector<Events> m_movePerIdx;

 public:
  enum Verbosity {
    TXT_NONE = 0,
    TXT_PER_DOMAIN = 1,
    TXT_PER_TYPE = 2,
    TXT_ALL = TXT_PER_DOMAIN | TXT_PER_TYPE
  };

  MoveStatRecorder() {}

  MoveStatRecorder(const MoveStatRecorder&) = delete;
  MoveStatRecorder operator=(const MoveStatRecorder&) = delete;

  void addEvent(const int idxDom, const std::string& domType,
                const bool moveSucceed, const bool moveAccepted) {
    if (int(m_movePerIdx.size()) <= idxDom) {
      m_movePerIdx.resize(idxDom + 1);
    }
    m_movePerIdx[idxDom].addEvent(moveSucceed, moveAccepted);
    m_movePerType[domType].addEvent(moveSucceed, moveAccepted);
  }

  std::string getString(const Verbosity inVerbosity) const {
    std::stringstream stream;

    stream << "[MOVE-STATS]\n";

    if (inVerbosity & TXT_PER_DOMAIN) {
      // Per type
      stream << "[MOVE-STATS] Per Domain:\n";
      for (size_t idx = 0; idx < m_movePerIdx.size(); ++idx) {
        stream << "[MOVE-STATS]   index = " << idx << "\n";
        stream << "[MOVE-STATS]     nb-moves = "
               << m_movePerIdx[idx].getNbTotalMoves() << "\n";
        stream << "[MOVE-STATS]     nb-succeed = "
               << m_movePerIdx[idx].getNbSucceedMoves() << " ("
               << m_movePerIdx[idx].getNbSucceedMovesPerc() << "%)\n";
        stream << "[MOVE-STATS]     nb-succeed-accepted = "
               << m_movePerIdx[idx].getNbAcceptedMoves() << " ("
               << m_movePerIdx[idx].getNbAcceptedMovesPerc() << "%)\n";
        stream << "[MOVE-STATS]     nb-failed = "
               << m_movePerIdx[idx].getNbFailMoves() << " ("
               << m_movePerIdx[idx].getNbFailMovesPerc() << "%)\n";
      }
    }

    // Per domain
    int totalNbSucceed = 0;
    int totalNbAccepted = 0;
    int totalNbFailed = 0;

    if (inVerbosity & TXT_PER_TYPE) {
      stream << "[MOVE-STATS] Per Type:\n";
    }
    for (const auto& iter : m_movePerType) {
      if (inVerbosity & TXT_PER_TYPE) {
        stream << "[MOVE-STATS]   type = " << iter.first << "\n";
        stream << "[MOVE-STATS]     nb-moves = "
               << iter.second.getNbTotalMoves() << "\n";
        stream << "[MOVE-STATS]     nb-succeed = "
               << iter.second.getNbSucceedMoves() << " ("
               << iter.second.getNbSucceedMovesPerc() << "%)\n";
        stream << "[MOVE-STATS]     nb-succeed-accepted = "
               << iter.second.getNbAcceptedMoves() << " ("
               << iter.second.getNbAcceptedMovesPerc() << "%)\n";
        stream << "[MOVE-STATS]     nb-failed = "
               << iter.second.getNbFailMoves() << " ("
               << iter.second.getNbFailMovesPerc() << "%)\n";
      }
      totalNbSucceed += iter.second.getNbSucceedMoves();
      totalNbAccepted += iter.second.getNbAcceptedMoves();
      totalNbFailed += iter.second.getNbFailMoves();
    }

    // For the simulation
    stream << "[MOVE-STATS] Total:\n";
    stream << "[MOVE-STATS]     nb-moves = " << totalNbSucceed + totalNbFailed
           << "\n";
    stream << "[MOVE-STATS]     nb-succeed = " << totalNbSucceed << " ("
           << 100. * double(totalNbSucceed) /
                  double(totalNbSucceed + totalNbFailed)
           << "%)\n";
    stream << "[MOVE-STATS]     nb-succeed-accepted = " << totalNbAccepted
           << " ("
           << 100. * double(totalNbAccepted) /
                  double(totalNbSucceed + totalNbFailed)
           << "%)\n";
    stream << "[MOVE-STATS]     nb-failed = " << totalNbFailed << " ("
           << 100. * double(totalNbFailed) /
                  double(totalNbSucceed + totalNbFailed)
           << "%)\n";

    return stream.str();
  }
};
}

#endif
