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
#ifndef EXCHANGELOGGER_H
#define EXCHANGELOGGER_H

#include <fmt/format.h>
#include <string>
#include <vector>

#include "util/log.h"
#include "util/util.h"

namespace mc {

class ExchangeLogger {
  struct Attempt {
    Attempt()
        : m_idx1(-1), m_idx2(-1), m_deltaEnergy(0), m_probability(0),
          m_exchanges(0), m_attempts(0) {}

    int m_idx1;
    int m_idx2;
    double m_deltaEnergy;
    double m_probability;
    int m_exchanges;
    int m_attempts;
  };

  int getAttemptStatsIdx(const int idx1, const int idx2) const {
    return std::max(idx1, idx2) * m_nbReplicas + std::min(idx1, idx2);
  }

  void print(const std::string &inStr, const int idxRep) {
    if (m_allLoggers[idxRep]) {
      DEBUG_ASSERT(m_allLoggers[idxRep] != nullptr,
                   "The log filename for replica {} has not been set.", idxRep);
      (*m_allLoggers[idxRep])(inStr);
    }
  }

  void printToAll(const std::string &inStr) {
    for (int idxRep = 0; idxRep < m_nbReplicas; ++idxRep) {
      print(inStr, idxRep);
    }
  }

  const int m_nbReplicas;
  std::vector<Attempt> m_attempts;
  std::vector<int> m_replicasAttemptIdx;
  int m_currentSweep;
  std::vector<Attempt> m_attemptsStats;
  int m_totalNbAttempts;
  std::vector<util::Logger *> m_allLoggers;

public:
  ExchangeLogger(const int inNbReplica)
      : m_nbReplicas(inNbReplica), m_currentSweep(-1), m_totalNbAttempts(0) {
    m_attemptsStats.resize(m_nbReplicas * m_nbReplicas);
    m_allLoggers.resize(m_nbReplicas, nullptr);
  }

  void setLoggerRef(const int idxRep, util::Logger &inLogger) {
    DEBUG_ASSERT(idxRep < m_nbReplicas,
                 "Replica idx must be lower than the number of replicas");
    m_allLoggers[idxRep] = &inLogger;
  }

  void printHeader(const std::vector<std::string> &configDirNames,
                   const std::string &configFilename, const int exchangeRate,
                   const int statisticRate, const int totalSweep) {
    DEBUG_ASSERT(
        m_nbReplicas == static_cast<int>(configDirNames.size()),
        "number of replicas from constructor and printHeader missmatch");

    std::string dirNamesStr;
    for (size_t idx = 0; idx < configDirNames.size(); ++idx) {
      dirNamesStr += fmt::format("id = {} : ", idx);
      dirNamesStr += configDirNames[idx];
      if (idx != configDirNames.size() - 1) {
        dirNamesStr += ", ";
      }
    }

    const std::string header =
        fmt::format("[Replica_Conf]\n"
                    "  number_replica: {}\n"
                    "  replica_directories: {}\n"
                    "  config_file_regex: {}\n"
                    "  exchange_rate: {}\n"
                    "  statistics_rate: {}\n"
                    "  total_mc_sweep: {}\n",
                    configDirNames.size(), dirNamesStr, configFilename,
                    exchangeRate, statisticRate, totalSweep);
    printToAll(header);
  }

  void startExchange(const int inCurrentSweep) {
    m_currentSweep = inCurrentSweep;
    m_attempts.clear();
    m_replicasAttemptIdx.clear();
    m_replicasAttemptIdx.resize(m_nbReplicas, -1);
  }

  void addAttempt(const int idx1, const int idx2, const double energy1,
                  const double energy2, const double probability,
                  const bool succeed) {
    DEBUG_ASSERT(m_replicasAttemptIdx[idx1] == -1,
                 "An attempt to exchange a replica must be done only once "
                 "(replica {}, sweep {})\n",
                 idx1, m_currentSweep);
    DEBUG_ASSERT(m_replicasAttemptIdx[idx2] == -1,
                 "An attempt to exchange a replica must be done only once "
                 "(replica {}, sweep {})\n",
                 idx2, m_currentSweep);

    m_totalNbAttempts += 1;

    Attempt currentAttempt;
    currentAttempt.m_idx1 = idx1;
    currentAttempt.m_idx2 = idx2;
    currentAttempt.m_deltaEnergy = (energy1 - energy2);
    currentAttempt.m_probability = probability;
    currentAttempt.m_exchanges = succeed ? 1 : 0;
    m_attempts.push_back(currentAttempt);

    const int idxCurrentAttempt = static_cast<int>(m_attempts.size() - 1);
    m_replicasAttemptIdx[idx1] = idxCurrentAttempt;
    m_replicasAttemptIdx[idx2] = idxCurrentAttempt;

    Attempt &currentAttemptStat =
        m_attemptsStats[getAttemptStatsIdx(idx1, idx2)];
    currentAttemptStat.m_probability =
        ((currentAttemptStat.m_probability *
          static_cast<double>(currentAttemptStat.m_attempts)) +
         probability) /
        static_cast<double>(currentAttemptStat.m_attempts + 1);
    currentAttemptStat.m_exchanges += succeed ? 1 : 0;
    currentAttemptStat.m_attempts += 1;
  }

  void printAttempts(const bool isOddEven) {
    // Exchanges str
    std::string exchangesStr;
    std::string exchangesProbabilities;
    if (isOddEven) {
      std::vector<std::pair<std::string, std::string>> neighborStats(
          m_nbReplicas - 1);

      for (size_t idxAttempt = 0; idxAttempt < m_attempts.size();
           ++idxAttempt) {
        const Attempt &attempt = m_attempts[idxAttempt];
        DEBUG_ASSERT(attempt.m_idx1 == attempt.m_idx2 - 1,
                     "Invalid attempts, must odd/even with idx1==idx2-1");
        neighborStats[attempt.m_idx1] = std::pair<std::string, std::string>(
            attempt.m_exchanges == 1 ? "x" : "o",
            fmt::format("{}", attempt.m_probability));
      }

      for (int idxRep = 0; idxRep < m_nbReplicas - 1; ++idxRep) {
        exchangesStr += fmt::format("{}-{}={}", idxRep, idxRep + 1,
                                    neighborStats[idxRep].first);
        exchangesProbabilities += fmt::format("{}-{}={}", idxRep, idxRep + 1,
                                              neighborStats[idxRep].second);
        if (idxRep != m_nbReplicas - 2) {
          exchangesStr += ", ";
          exchangesProbabilities += ", ";
        }
      }
    } else {
      for (size_t idxAttempt = 0; idxAttempt < m_attempts.size();
           ++idxAttempt) {
        const Attempt &attempt = m_attempts[idxAttempt];
        exchangesStr += fmt::format("{}-{}={}", attempt.m_idx1, attempt.m_idx2,
                                    attempt.m_exchanges == 1 ? "x" : "o");
        exchangesProbabilities += fmt::format(
            "{}-{}={}", attempt.m_idx1, attempt.m_idx2, attempt.m_probability);
        if (idxAttempt != m_attempts.size() - 1) {
          exchangesStr += ", ";
          exchangesProbabilities += ", ";
        }
      }
    }

    for (int idxRep = 0; idxRep < m_nbReplicas; ++idxRep) {
      std::string repDeltaEnergy;
      std::string repAttmptIdxs;
      if (m_replicasAttemptIdx[idxRep] != -1) {
        repDeltaEnergy =
            fmt::format("  delta_energy: {} (kT)\n",
                        m_attempts[m_replicasAttemptIdx[idxRep]].m_deltaEnergy);
        repAttmptIdxs =
            fmt::format("  attempted_exchange: {}-{}\n",
                        m_attempts[m_replicasAttemptIdx[idxRep]].m_idx1,
                        m_attempts[m_replicasAttemptIdx[idxRep]].m_idx2);
      }

      const std::string toprint =
          fmt::format("[Replica_Attempt]\n"
                      "  sweep: {}\n"
                      "{}"
                      "{}"
                      "  exchanges: {}\n"
                      "  probabilities: {}\n",
                      m_currentSweep, repDeltaEnergy, repAttmptIdxs,
                      exchangesStr, exchangesProbabilities);

      print(toprint, idxRep);
    }
  }

  void printStats() {
    std::string exchangesProbabilities;
    std::string exchangesNumber;
    std::string exchangesAverage;

    std::vector<double> transitionsMatrix(m_nbReplicas * m_nbReplicas, 0);
    for (int idxAttempt = 0; idxAttempt < m_nbReplicas; ++idxAttempt) {
      transitionsMatrix[getAttemptStatsIdx(idxAttempt, idxAttempt)] = 1.;
    }

    const std::string separator = ", ";

    for (int idx1 = 0; idx1 < m_nbReplicas; ++idx1) {
      for (int idx2 = idx1 + 1; idx2 < m_nbReplicas; ++idx2) {
        const Attempt &currentAttemptStat =
            m_attemptsStats[getAttemptStatsIdx(idx1, idx2)];
        if (currentAttemptStat.m_attempts) {
          exchangesProbabilities +=
              fmt::format("{}-{}={}{}", idx1, idx2,
                          currentAttemptStat.m_probability, separator);
          exchangesNumber +=
              fmt::format("{}-{}={}{}", idx1, idx2,
                          currentAttemptStat.m_exchanges, separator);
          exchangesAverage += fmt::format(
              "{}-{}={}{}", idx1, idx2,
              static_cast<double>(currentAttemptStat.m_exchanges) /
                  static_cast<double>(currentAttemptStat.m_attempts),
              separator);

          if (currentAttemptStat.m_exchanges) {
            const double transitionCoef =
                static_cast<double>(currentAttemptStat.m_exchanges) /
                static_cast<double>(m_totalNbAttempts);
            transitionsMatrix[getAttemptStatsIdx(idx1, idx1)] -= transitionCoef;
            transitionsMatrix[getAttemptStatsIdx(idx2, idx2)] -= transitionCoef;
            transitionsMatrix[getAttemptStatsIdx(idx1, idx2)] += transitionCoef;
          }
        }
      }
    }

    if (m_totalNbAttempts) {
      exchangesProbabilities.resize(exchangesProbabilities.size() -
                                    separator.size());
      exchangesNumber.resize(exchangesNumber.size() - separator.size());
      exchangesAverage.resize(exchangesAverage.size() - separator.size());
    }

    std::string stats = fmt::format("[Replica_Average]\n"
                                    "  total_attempts: {}\n"
                                    "  average_probabilities: {}\n"
                                    "  number_exchange: {}\n"
                                    "  average_exchange: {}\n"
                                    "  empirical_transition_matrix:\n",
                                    m_totalNbAttempts, exchangesProbabilities,
                                    exchangesNumber, exchangesAverage);

    for (int idx1 = 0; idx1 < m_nbReplicas; ++idx1) {
      stats += "   ";
      for (int idx2 = 0; idx2 < m_nbReplicas; ++idx2) {
        stats += fmt::format("{:06.5f} ",
                             transitionsMatrix[getAttemptStatsIdx(idx1, idx2)]);
      }
      stats += "\n";
    }

    printToAll(stats);
  }
};
} // namespace mc

#endif
