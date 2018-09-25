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
#ifndef ENERGY_ENERGYMATRIX_H
#define ENERGY_ENERGYMATRIX_H

#include <array>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <iostream>
#include <omp.h>
#include <vector>

#include "io/serializer.h"
#include "util/util.h"

namespace energy {

template <typename T, int BlockingSize = 32>
class EnergyMatrix : public io::AbstractSerializable {
  static_assert(BlockingSize > 0, "BlockingSize must geater than 0");

  static constexpr int NbMainValues =
      3; // for energy + inner energy + energy for connections
  static constexpr int InnerEnergyIdx = 1;
  static constexpr int ConnectionEnergyIdx = 2;

  int m_nbRowDomains;
  int m_nbColDomains;
  int m_nbContributions;
  int m_nbValues;
  int m_nbRowBlocks;
  int m_nbColBlocks;
  std::vector<T> m_pData;

  int indexOf(const int inIdxRow, const int inIdxCol) const {
    const int idxRowBlock = inIdxRow / BlockingSize;
    const int idxColBlock = inIdxCol / BlockingSize;
    const int idxValuesBlock = (idxRowBlock + idxColBlock * m_nbRowBlocks) *
                               BlockingSize * BlockingSize;
    return idxValuesBlock + (inIdxRow % BlockingSize) +
           (inIdxCol % BlockingSize) * BlockingSize;
  }

  T &getEnergyWrite(const int idDomainRow, const int idDomainCol) {
    DEBUG_ASSERT(idDomainRow < m_nbRowDomains && idDomainCol < m_nbColDomains,
                 "Out of bound array access ({} {}) with rows={}, cols={}",
                 idDomainRow, idDomainCol, m_nbRowDomains, m_nbColDomains);
    return m_pData[indexOf(idDomainRow, idDomainCol) * m_nbValues];
  }

  T &getEnergyConnectionsWrite(const int idDomainRow, const int idDomainCol) {
    DEBUG_ASSERT(idDomainRow < m_nbRowDomains && idDomainCol < m_nbColDomains,
                 "Out of bound array access ({} {}) with rows={}, cols={}",
                 idDomainRow, idDomainCol, m_nbRowDomains, m_nbColDomains);
    return m_pData[indexOf(idDomainRow, idDomainCol) * m_nbValues +
                   ConnectionEnergyIdx];
  }

  T &getContributionWrite(const int idDomainRow, const int idDomainCol,
                          const int idxContribution) {
    DEBUG_ASSERT(idDomainRow < m_nbRowDomains && idDomainCol < m_nbColDomains &&
                     idxContribution < m_nbContributions,
                 "Out of bound array access ({} {} {}) with rows={}, cols={}, "
                 "contributions={}",
                 idDomainRow, idDomainCol, idxContribution, m_nbRowDomains,
                 m_nbColDomains, m_nbContributions);
    return m_pData[indexOf(idDomainRow, idDomainCol) * m_nbValues +
                   NbMainValues + idxContribution];
  }

public:
  static const bool IsRowMajor = true;

  explicit EnergyMatrix(int inNbRowDomains_, int inNbColDomains_,
                        int nbContributions_)
      : m_nbRowDomains(inNbRowDomains_), m_nbColDomains(inNbColDomains_),
        m_nbContributions(nbContributions_),
        m_nbValues(nbContributions_ + NbMainValues),
        m_nbRowBlocks((m_nbRowDomains + BlockingSize - 1) / BlockingSize),
        m_nbColBlocks((m_nbColDomains + BlockingSize - 1) / BlockingSize) {
    m_pData.resize(m_nbRowBlocks * m_nbColBlocks * BlockingSize * BlockingSize *
                   m_nbValues);
  }
  explicit EnergyMatrix(std::size_t inNbRowDomains_,
                        std::size_t inNbColDomains_, int nbValues_)
      : EnergyMatrix(static_cast<int>(inNbRowDomains_),
                     static_cast<int>(inNbColDomains_), nbValues_) {}
  // support moving
  EnergyMatrix(EnergyMatrix &&rhs) = default;
  EnergyMatrix &operator=(EnergyMatrix &&rhs) = default;
  // support copying
  EnergyMatrix(const EnergyMatrix &rhs) = default;
  EnergyMatrix &operator=(const EnergyMatrix &rhs) = default;

  void serialize(io::Serializer &serializer) const final {
    serializer.append(m_nbRowDomains, "m_nbRowDomains");
    serializer.append(m_nbColDomains, "m_nbColDomains");
    serializer.append(m_nbContributions, "m_nbContributions");
    serializer.append(m_nbValues, "m_nbValues");
    serializer.append(m_nbRowBlocks, "m_nbRowBlocks");
    serializer.append(m_nbColBlocks, "m_nbColBlocks");
    serializer.append(m_pData, "m_pData");
  }

  EnergyMatrix(io::Deserializer &deserializer)
      : m_nbRowDomains(
            deserializer.restore<decltype(m_nbRowDomains)>("m_nbRowDomains")),
        m_nbColDomains(
            deserializer.restore<decltype(m_nbColDomains)>("m_nbColDomains")),
        m_nbContributions(deserializer.restore<decltype(m_nbContributions)>(
            "m_nbContributions")),
        m_nbValues(deserializer.restore<decltype(m_nbValues)>("m_nbValues")),
        m_nbRowBlocks(
            deserializer.restore<decltype(m_nbRowBlocks)>("m_nbRowBlocks")),
        m_nbColBlocks(
            deserializer.restore<decltype(m_nbColBlocks)>("m_nbColBlocks")),
        m_pData(deserializer.restore<decltype(m_pData)>("m_pData")) {
    DEBUG_ASSERT(static_cast<int>(m_pData.size()) ==
                     m_nbRowBlocks * m_nbColBlocks * BlockingSize *
                         BlockingSize * m_nbValues,
                 "Invalid size");
  }

  /////////////////////////////////////////////////////////////////////////////

  const T &getEnergy(const int idDomainRow, const int idDomainCol) const {
    DEBUG_ASSERT(idDomainRow < m_nbRowDomains && idDomainCol < m_nbColDomains,
                 "Out of bound array access ({} {}) with rows={}, cols={}",
                 idDomainRow, idDomainCol, m_nbRowDomains, m_nbColDomains);
    return m_pData[indexOf(idDomainRow, idDomainCol) * m_nbValues];
  }

  /////////////////////////////////////////////////////////////////////////////

  void addEnergyConnections(const int idDomainRow, const int idDomainCol,
                            const T &inVal) {
    DEBUG_ASSERT(idDomainRow < m_nbRowDomains && idDomainCol < m_nbColDomains,
                 "Out of bound array access ({} {}) with rows={}, cols={}",
                 idDomainRow, idDomainCol, m_nbRowDomains, m_nbColDomains);
    T &target = m_pData[indexOf(idDomainRow, idDomainCol) * m_nbValues +
                        ConnectionEnergyIdx];
    target += inVal;
    T &targetEnergy = m_pData[indexOf(idDomainRow, idDomainCol) * m_nbValues];
    targetEnergy += inVal;
  }

  const T &getEnergyConnections(const int idDomainRow,
                                const int idDomainCol) const {
    DEBUG_ASSERT(idDomainRow < m_nbRowDomains && idDomainCol < m_nbColDomains,
                 "Out of bound array access ({} {}) with rows={}, cols={}",
                 idDomainRow, idDomainCol, m_nbRowDomains, m_nbColDomains);
    return m_pData[indexOf(idDomainRow, idDomainCol) * m_nbValues +
                   ConnectionEnergyIdx];
  }

  /////////////////////////////////////////////////////////////////////////////

  void addContribution(const int idDomainRow, const int idDomainCol,
                       const int idxContribution, const T &inVal) {
    DEBUG_ASSERT(idDomainRow < m_nbRowDomains && idDomainCol < m_nbColDomains &&
                     idxContribution < m_nbContributions,
                 "Out of bound array access ({} {} {}) with rows={}, cols={}, "
                 "contributions={}",
                 idDomainRow, idDomainCol, idxContribution, m_nbRowDomains,
                 m_nbColDomains, m_nbContributions);
    T &target = m_pData[indexOf(idDomainRow, idDomainCol) * m_nbValues +
                        NbMainValues + idxContribution];
    target += inVal;
    T &targetEnergy = m_pData[indexOf(idDomainRow, idDomainCol) * m_nbValues];
    targetEnergy += inVal;
  }

  const T &getContribution(const int idDomainRow, const int idDomainCol,
                           const int idxContribution) const {
    DEBUG_ASSERT(idDomainRow < m_nbRowDomains && idDomainCol < m_nbColDomains &&
                     idxContribution < m_nbContributions,
                 "Out of bound array access ({} {} {}) with rows={}, cols={}, "
                 "contributions={}",
                 idDomainRow, idDomainCol, idxContribution, m_nbRowDomains,
                 m_nbColDomains, m_nbContributions);
    return m_pData[indexOf(idDomainRow, idDomainCol) * m_nbValues +
                   NbMainValues + idxContribution];
  }

  /////////////////////////////////////////////////////////////////////////////

  int nbRowDomains() const noexcept { return m_nbRowDomains; };
  int nbColDomains() const noexcept { return m_nbColDomains; };
  int nbEnergyValues() const noexcept { return m_nbValues; };
  int nbContributions() const noexcept { return m_nbContributions; };

  /////////////////////////////////////////////////////////////////////////////

  void reset() { setAll(T()); }

  void setAll(const T &inVal) {
    for (size_t idx = 0; idx < m_pData.size(); ++idx) {
      m_pData[idx] = inVal;
    }
  }

  bool operator==(const EnergyMatrix &other) const {
    return m_nbRowDomains == other.m_nbRowDomains &&
           m_nbColDomains == other.m_nbColDomains &&
           m_nbValues == other.m_nbValues &&
           m_nbContributions == other.m_nbContributions &&
           m_pData == other.m_pData;
  }

  bool operator!=(const EnergyMatrix &other) const { return !(*this == other); }

  bool sameDimension(const EnergyMatrix &other) const {
    return m_nbRowDomains == other.m_nbRowDomains &&
           m_nbColDomains == other.m_nbColDomains &&
           m_nbValues == other.m_nbValues &&
           m_nbContributions == other.m_nbContributions;
  }

  /////////////////////////////////////////////////////////////////////////////

  void replaceRowAndCol(const int rowColIdx, const EnergyMatrix &inValues) {
    replaceRow(rowColIdx, inValues);
    replaceCol(rowColIdx, inValues);
  }

  void replaceRow(const int rowIdx, const EnergyMatrix &inRowValues) {
    DEBUG_ASSERT(inRowValues.m_nbValues == m_nbValues &&
                     inRowValues.m_nbContributions == m_nbContributions &&
                     ((inRowValues.m_nbRowDomains == 1 &&
                       inRowValues.m_nbColDomains == m_nbColDomains) ||
                      (inRowValues.m_nbRowDomains == m_nbRowDomains &&
                       inRowValues.m_nbColDomains == 1)),
                 "Vector has invalid size {} {} {} {} (should be compatible "
                 "with {} {} {} {}",
                 inRowValues.m_nbRowDomains, inRowValues.m_nbColDomains,
                 inRowValues.m_nbValues, inRowValues.m_nbContributions,
                 m_nbRowDomains, m_nbColDomains, m_nbValues, m_nbContributions);
    if (inRowValues.m_nbColDomains == m_nbColDomains) {
      for (int idxCol = 0; idxCol < m_nbColDomains; ++idxCol) {
        (*this).getEnergyWrite(rowIdx, idxCol) =
            inRowValues.getEnergy(0, idxCol);
        (*this).getEnergyConnectionsWrite(rowIdx, idxCol) =
            inRowValues.getEnergyConnections(0, idxCol);
        for (int idxVal = 0; idxVal < m_nbContributions; ++idxVal) {
          (*this).getContributionWrite(rowIdx, idxCol, idxVal) =
              inRowValues.getContribution(0, idxCol, idxVal);
        }
      }
    } else {
      for (int idxCol = 0; idxCol < m_nbColDomains; ++idxCol) {
        (*this).getEnergyWrite(rowIdx, idxCol) =
            inRowValues.getEnergy(idxCol, 0);
        (*this).getEnergyConnectionsWrite(rowIdx, idxCol) =
            inRowValues.getEnergyConnections(idxCol, 0);
        for (int idxVal = 0; idxVal < m_nbContributions; ++idxVal) {
          (*this).getContributionWrite(rowIdx, idxCol, idxVal) =
              inRowValues.getContribution(idxCol, 0, idxVal);
        }
      }
    }
  }

  void replaceCol(const int colIdx, const EnergyMatrix &inColValues) {
    DEBUG_ASSERT(inColValues.m_nbValues == m_nbValues &&
                     inColValues.m_nbContributions == m_nbContributions &&
                     ((inColValues.m_nbRowDomains == 1 &&
                       inColValues.m_nbColDomains == m_nbColDomains) ||
                      (inColValues.m_nbRowDomains == m_nbRowDomains &&
                       inColValues.m_nbColDomains == 1)),
                 "Vector has invalid size {} {} {} {} (should be compatible "
                 "with {} {} {} {}",
                 inColValues.m_nbRowDomains, inColValues.m_nbColDomains,
                 inColValues.m_nbValues, inColValues.m_nbContributions,
                 m_nbRowDomains, m_nbColDomains, m_nbValues, m_nbContributions);
    if (inColValues.m_nbRowDomains == m_nbRowDomains) {
      for (int idxRow = 0; idxRow < m_nbRowDomains; ++idxRow) {
        (*this).getEnergyWrite(idxRow, colIdx) =
            inColValues.getEnergy(idxRow, 0);
        (*this).getEnergyConnectionsWrite(idxRow, colIdx) =
            inColValues.getEnergyConnections(idxRow, 0);
        for (int idxVal = 0; idxVal < m_nbContributions; ++idxVal) {
          (*this).getContributionWrite(idxRow, colIdx, idxVal) =
              inColValues.getContribution(idxRow, 0, idxVal);
        }
      }
    } else {
      for (int idxRow = 0; idxRow < m_nbRowDomains; ++idxRow) {
        (*this).getEnergyWrite(idxRow, colIdx) =
            inColValues.getEnergy(0, idxRow);
        (*this).getEnergyConnectionsWrite(idxRow, colIdx) =
            inColValues.getEnergyConnections(0, idxRow);
        for (int idxVal = 0; idxVal < m_nbContributions; ++idxVal) {
          (*this).getContributionWrite(idxRow, colIdx, idxVal) =
              inColValues.getContribution(0, idxRow, idxVal);
        }
      }
    }
  }

  /////////////////////////////////////////////////////////////////////////////

  T getTotalEnergy() const {
    T sum = 0;
    for (auto j = 0; j < nbColDomains(); ++j) {
      for (auto i = j; i < nbRowDomains(); ++i) {
        sum += getEnergy(i, j);
      }
    }
#ifndef NDEBUG
    checkCoherency();
#endif
    return sum;
  }

  T getTotalEnergyConnections() const {
    T sum = 0;
    for (auto j = 0; j < nbColDomains(); ++j) {
      for (auto i = j; i < nbRowDomains(); ++i) {
        sum += getEnergyConnections(i, j);
      }
    }
    return sum;
  }

  std::vector<T> getTotalContributions() const {
    std::vector<T> sum(nbContributions());
    for (auto j = 0; j < nbColDomains(); ++j) {
      for (auto i = j; i < nbRowDomains(); ++i) {
        for (auto v = 0; v < nbContributions(); ++v) {
          sum[v] += getContribution(i, j, v);
        }
      }
    }
    return sum;
  }

  void checkCoherency() const {
    for (auto j = 0; j < nbColDomains(); ++j) {
      for (auto i = j; i < nbRowDomains(); ++i) {
        double sumContributions = 0;
        for (auto v = 0; v < nbContributions(); ++v) {
          sumContributions += getContribution(i, j, v);
        }
        const double energyConnections = getEnergyConnections(i, j);
        const double sumComponents = sumContributions + energyConnections;

        const double totalEnergy = getEnergy(i, j);

        const double diff =
            (totalEnergy == 0 ? sumComponents
                              : (std::abs(sumComponents - totalEnergy) /
                                 (std::abs(totalEnergy) +
                                  std::numeric_limits<double>::epsilon())));
        if (diff >= 1e-8) {
          throw std::runtime_error(fmt::format(
              "EnergyMatrix -- Error the contributions are not "
              "equal to the energy, should be {}, is {} (diff = "
              "{}, connections = {}, contributions = {}), for indexes {} {}",
              totalEnergy, sumComponents, diff, energyConnections,
              sumContributions, i, j));
        }
      }
    }
  }
};

using rEnergyMatrix = EnergyMatrix<double>;

template <int BufferSize = 1024> class EnergyMatrixBuffer {
  static_assert(BufferSize > 0, "BufferSize must be positive");

  struct BufferValue {
    int idxDomRow;
    int idxDomCol;
    int inOffset;
    int inIdx;
    double inVal;
  };

  rEnergyMatrix &m_energyArray;
  BufferValue m_values[BufferSize];
  int m_nbValues;

  omp_lock_t &m_energyMatrixMutex;

public:
  EnergyMatrixBuffer(rEnergyMatrix &inEnergyArray,
                     omp_lock_t &inEnergyMatrixMutex)
      : m_energyArray(inEnergyArray), m_nbValues(0),
        m_energyMatrixMutex(inEnergyMatrixMutex) {}

  ~EnergyMatrixBuffer() { flushWithLock(); }

  EnergyMatrixBuffer(const EnergyMatrixBuffer &) = delete;
  EnergyMatrixBuffer(EnergyMatrixBuffer &&) = delete;
  EnergyMatrixBuffer &operator=(const EnergyMatrixBuffer &) = delete;
  EnergyMatrixBuffer &operator=(EnergyMatrixBuffer &&) = delete;

  EnergyMatrixBuffer &addContribution(const int idxDomRow, const int idxDomCol,
                                      const int inOffset, const int inIdx,
                                      const double inVal) {
    if (m_nbValues == BufferSize) {
      flushWithLock();
    }
    DEBUG_ASSERT(inOffset + inIdx <
                     static_cast<int>(m_energyArray.nbEnergyValues()),
                 "Invalid inIdx {} + inOffset {} >= m_contributions.size() {}",
                 inIdx, inOffset, m_energyArray.nbEnergyValues());
    m_values[m_nbValues++] =
        BufferValue{idxDomRow, idxDomCol, inOffset, inIdx, inVal};
    return *this;
  }

  void flushWithLock() {
    omp_set_lock(&m_energyMatrixMutex);
    for (int idxValue = 0; idxValue < m_nbValues; ++idxValue) {
      m_energyArray.addContribution(
          m_values[idxValue].idxDomRow, m_values[idxValue].idxDomCol,
          m_values[idxValue].inOffset + m_values[idxValue].inIdx,
          m_values[idxValue].inVal);
    }
    m_nbValues = 0;
    omp_unset_lock(&m_energyMatrixMutex);
  }

  class Accesser {
    EnergyMatrixBuffer &m_buffer;
    const int m_idxDomRow;
    const int m_idxDomCol;

    Accesser(EnergyMatrixBuffer &inBuffer, const int idxDomRow,
             const int idxDomCol)
        : m_buffer(inBuffer), m_idxDomRow(idxDomRow), m_idxDomCol(idxDomCol) {}

  public:
    Accesser(const Accesser &) = default;
    Accesser(Accesser &&) = default;
    Accesser &operator=(const Accesser &) = default;
    Accesser &operator=(Accesser &&) = default;

    Accesser &addContribution(const int inOffset, const int inIdx,
                              const double inVal) {
      m_buffer.addContribution(m_idxDomRow, m_idxDomCol, inOffset, inIdx,
                               inVal);
      return *this;
    }

    friend EnergyMatrixBuffer;
  };

  Accesser getAccesser(const int idxDomRow, const int idxDomCol) {
    return Accesser(*this, idxDomRow, idxDomCol);
  }
};

} // namespace energy

namespace std {
template <typename T>
std::ostream &operator<<(std::ostream &out,
                         const energy::EnergyMatrix<T> &mat) {
  const auto nbDomains = mat.nbRowDomains();
  out << fmt::format("Energy Mat:\n");
  for (auto i = 0; i < nbDomains; ++i) {
    out << fmt::format("\t");
    for (auto j = 0; j < nbDomains - 1; ++j) {
      out << fmt::format("{:12.8f}, ", mat.getEnergy(i, j));
    }
    out << fmt::format("{:12.8f}\n", mat.getEnergy(i, nbDomains - 1));
  }
  return out;
}
} // namespace std
#endif // ENERGY_ENERGYMATRIX_H
