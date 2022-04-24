#ifndef PAIRKERNELMANAGER_H
#define PAIRKERNELMANAGER_H

#include <memory>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

#include "abstractpairkernel.h"
#include "domains/abstractdomain.h"
#include "energy/forcefield.h"
#include "pairkernelbuilder.h"
#include "util/array.h"
#include "util/util.h"

namespace pairkernels {

class PairKernelManager : public io::AbstractSerializable {
  std::vector<std::unique_ptr<AbstractPairKernel>> m_kernels;
  util::Array<AbstractPairKernel*> m_kernelsMapping;

  int m_nbContributions;

 public:
  void serialize(io::Serializer& serializer) const final {
    serializer.append(m_nbContributions, "m_nbContributions");

    serializer.append(m_kernelsMapping.rows(), "m_kernelsMapping.rows");
    serializer.append(m_kernelsMapping.cols(), "m_kernelsMapping.cols");

    int checksum = 0;

    for (int idxCol = 0; idxCol < m_kernelsMapping.cols(); ++idxCol) {
      for (int idxRow = 0; idxRow < m_kernelsMapping.rows(); ++idxRow) {
        if (m_kernelsMapping(idxRow, idxCol)) {
          serializer.append(int(1), "next");
          serializer.append(idxRow, "idxRow");
          serializer.append(idxCol, "idxCol");
          serializer.append(m_kernelsMapping(idxRow, idxCol)->type(), "type");
          checksum += 1;
        }
      }
    }
    serializer.append(int(0), "next");
    serializer.append(checksum, "checksum");
  }

  PairKernelManager(io::Deserializer& deserializer)
      : m_kernels(0),
        m_kernelsMapping(0, 0),
        m_nbContributions(deserializer.restore<decltype(m_nbContributions)>(
            "m_nbContributions")) {
    const int m_kernelsMapping_rows =
        deserializer.restore<decltype(m_kernelsMapping.rows())>(
            "m_kernelsMapping.rows");
    const int m_kernelsMapping_cols =
        deserializer.restore<decltype(m_kernelsMapping.cols())>(
            "m_kernelsMapping.cols");

    m_kernelsMapping = util::Array<AbstractPairKernel*>(
        m_kernelsMapping_rows, m_kernelsMapping_cols, nullptr);

    int checksum = 0;

    while (deserializer.restore<int>("next")) {
      const int idxRow = deserializer.restore<int>("idxRow");
      const int idxCol = deserializer.restore<int>("idxCol");
      const std::string type = deserializer.restore<std::string>("type");
      const auto& kernel =
          std::find_if(std::begin(m_kernels), std::end(m_kernels),
                       [&type](const std::unique_ptr<AbstractPairKernel>& k) {
                         return k->type() == type;
                       });
      if (kernel == std::end(m_kernels)) {
        m_kernels.emplace_back(PairKernelBuilder::Build(type));
        m_kernelsMapping(idxRow, idxCol) = m_kernels.back().get();
      } else {
        m_kernelsMapping(idxRow, idxCol) = kernel->get();
      }
      checksum += 1;
    }
    const int recvCheckSum = deserializer.restore<int>("checksum");

    if (checksum != recvCheckSum) {
      throw std::runtime_error(
          fmt::format("Invalid number of kernel unpack {}, should be {}",
                      checksum, recvCheckSum));
    }
  }

  explicit PairKernelManager()
      : m_kernels(0), m_kernelsMapping(0, 0), m_nbContributions(0) {}

  explicit PairKernelManager(
      const std::vector<std::array<std::string, 3>>& inAllPairs,
      const domains::Domains& topology)
      : m_kernels(0), m_kernelsMapping(0, 0), m_nbContributions(0) {
    const std::string defaultKey = "default";
    const std::string allKey = "*";

    std::unordered_map<std::string, std::unique_ptr<AbstractPairKernel>>
        kernels;
    std::map<std::string, int> allTypesIds;
    std::map<int, std::string> allTypesIdsStr;
    int maxTypeId = 0;
    for (auto& dom : topology) {
      if (allTypesIds.find(dom->name()) != allTypesIds.end() &&
          allTypesIds[dom->name()] != dom->typeId()) {
        throw std::invalid_argument(
            "Some domains are from the same type but have different type ID");
      }
      if (allTypesIdsStr.find(dom->typeId()) != allTypesIdsStr.end() &&
          allTypesIdsStr[dom->typeId()] != dom->name()) {
        throw std::invalid_argument(
            "Some domains are from the different type but have the same type "
            "ID");
      }
      allTypesIds[dom->name()] = dom->typeId();
      allTypesIdsStr[dom->typeId()] = dom->name();
      maxTypeId = std::max(maxTypeId, dom->typeId());
    }

    m_kernelsMapping =
        util::Array<AbstractPairKernel*>(maxTypeId + 1, maxTypeId + 1, nullptr);

    util::Array<int> defaultCanBeRemoved(maxTypeId + 1, maxTypeId + 1, true);

    for (auto& pair : inAllPairs) {
      // if both pairs are */Default
      if ((pair[0] == allKey || pair[0] == defaultKey) &&
          (pair[1] == allKey || pair[1] == defaultKey)) {
        const bool setByDefault =
            (pair[0] == defaultKey && pair[1] == defaultKey);
        if (kernels.find(pair[2]) == kernels.end()) {
          kernels[pair[2]] = PairKernelBuilder::Build(pair[2]);
        }
        for (auto other1 : allTypesIds) {
          for (auto other2 : allTypesIds) {
            if (defaultCanBeRemoved(other1.second, other2.second) == false &&
                setByDefault == false &&
                kernelExists(other1.second, other2.second)) {
              throw std::invalid_argument(fmt::format(
                  "PairKernelManager -- pair {} / {} has already been defined",
                  other1.first, other2.first));
            }
            if (setByDefault == false ||
                defaultCanBeRemoved(other1.second, other2.second)) {
              m_kernelsMapping(other1.second, other2.second) =
                  kernels[pair[2]].get();
              defaultCanBeRemoved(other1.second, other2.second) = setByDefault;
            }
          }
        }
      }
      // If one is */Default and the other exist
      else if (((pair[0] == allKey || pair[0] == defaultKey) &&
                allTypesIds.find(pair[1]) != allTypesIds.end()) ||
               ((pair[1] == allKey || pair[1] == defaultKey) &&
                allTypesIds.find(pair[0]) != allTypesIds.end())) {
        const int globalIdx =
            ((pair[0] == allKey || pair[0] == defaultKey) ? 0 : 1);
        const std::string& globalStr = pair[globalIdx];
        const std::string& specificType = pair[!globalIdx];
        const int specificTypeId = allTypesIds[specificType];

        const bool setByDefault = (globalStr == defaultKey);
        if (kernels.find(pair[2]) == kernels.end()) {
          kernels[pair[2]] = PairKernelBuilder::Build(pair[2]);
        }
        for (auto other : allTypesIds) {
          if ((defaultCanBeRemoved(other.second, specificTypeId) == false ||
               defaultCanBeRemoved(specificTypeId, other.second) == false) &&
              setByDefault == false &&
              kernelExists(other.second, specificTypeId)) {
            throw std::invalid_argument(fmt::format(
                "PairKernelManager -- pair {} / {} has already been defined",
                other.first, specificType));
          }
          if (setByDefault == false ||
              defaultCanBeRemoved(other.second, specificTypeId)) {
            m_kernelsMapping(other.second, specificTypeId) =
                kernels[pair[2]].get();
            m_kernelsMapping(specificTypeId, other.second) =
                kernels[pair[2]].get();
            defaultCanBeRemoved(other.second, specificTypeId) = setByDefault;
            defaultCanBeRemoved(specificTypeId, other.second) = setByDefault;
          }
        }
      } else if (pair[0] != allKey && pair[0] != defaultKey &&
                 allTypesIds.find(pair[0]) != allTypesIds.end() &&
                 pair[1] != allKey && pair[1] != defaultKey &&
                 allTypesIds.find(pair[1]) != allTypesIds.end()) {
        if ((defaultCanBeRemoved(allTypesIds[pair[0]], allTypesIds[pair[1]]) ==
                 false ||
             defaultCanBeRemoved(allTypesIds[pair[1]], allTypesIds[pair[0]]) ==
                 false) &&
            kernelExists(allTypesIds[pair[0]], allTypesIds[pair[1]])) {
          throw std::invalid_argument(fmt::format(
              "PairKernelManager -- pair {} / {} has already been defined",
              pair[0], pair[1]));
        }
        if (kernels.find(pair[2]) == kernels.end()) {
          kernels[pair[2]] = PairKernelBuilder::Build(pair[2]);
        }
        m_kernelsMapping(allTypesIds[pair[0]], allTypesIds[pair[1]]) =
            kernels[pair[2]].get();
        m_kernelsMapping(allTypesIds[pair[1]], allTypesIds[pair[0]]) =
            kernels[pair[2]].get();
        defaultCanBeRemoved(allTypesIds[pair[0]], allTypesIds[pair[1]]) = false;
        defaultCanBeRemoved(allTypesIds[pair[1]], allTypesIds[pair[0]]) = false;
      }
    }

    std::set<const AbstractPairKernel*> usedKernels;
    for (auto other1 : allTypesIds) {
      for (auto other2 : allTypesIds) {
        if (!kernelExists(other1.second, other2.second)) {
          throw std::invalid_argument(fmt::format(
              "PairKernelManager -- pair {} / {} has never been defined",
              other1.first, other2.first));
        }
        usedKernels.insert(&getKernel(other1.second, other2.second));
      }
    }

    // Remove unused kernels from the map
    const auto kernelIterEnd = kernels.end();
    auto kernelIter = kernels.begin();
    while (kernelIter != kernelIterEnd) {
      if (usedKernels.find(kernelIter->second.get()) == usedKernels.end()) {
        kernelIter = kernels.erase(kernelIter);
      } else {
        ++kernelIter;
      }
    }

    m_nbContributions = 0;  // We do not count energy
    for (auto& kernel : kernels) {
      kernel.second->setOffsetContributions(m_nbContributions);
      m_nbContributions += kernel.second->getNbExtraContributions();
      m_kernels.push_back(std::move(kernel.second));
    }
  }

  PairKernelManager(const PairKernelManager&) = delete;
  PairKernelManager& operator=(const PairKernelManager&) = delete;

  PairKernelManager(PairKernelManager&& inOther)
      : m_kernels(std::move(inOther.m_kernels)),
        m_kernelsMapping(std::move(inOther.m_kernelsMapping)),
        m_nbContributions(std::move(inOther.m_nbContributions)) {}

  PairKernelManager& operator=(PairKernelManager&& inOther) {
    m_nbContributions = std::move(inOther.m_nbContributions);
    m_kernelsMapping = std::move(inOther.m_kernelsMapping);
    m_kernels = std::move(inOther.m_kernels);
    return *this;
  }

  bool kernelExists(const int inDomTypeId1, const int inDomTypeId2) const {
    DEBUG_ASSERT(inDomTypeId1 < m_kernelsMapping.rows(), "Type id too large");
    DEBUG_ASSERT(inDomTypeId2 < m_kernelsMapping.cols(), "Type id too large");
    return m_kernelsMapping(inDomTypeId1, inDomTypeId2) != nullptr;
  }

  const AbstractPairKernel& getKernel(const int inDomTypeId1,
                                      const int inDomTypeId2) const {
    DEBUG_ASSERT(inDomTypeId1 < m_kernelsMapping.rows(), "Type id too large");
    DEBUG_ASSERT(inDomTypeId2 < m_kernelsMapping.cols(), "Type id too large");
    DEBUG_ASSERT(m_kernelsMapping(inDomTypeId1, inDomTypeId2) != nullptr,
                 "Kernel do not exist for ids {} , {}", inDomTypeId1,
                 inDomTypeId2);
    return *m_kernelsMapping(inDomTypeId1, inDomTypeId2);
  }

  //////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////

  int getNbContributions() const { return m_nbContributions; }

  std::vector<std::string> getContributionLabels() const {
    std::vector<std::string> labels;
    for (const auto& kernel : m_kernels) {
      const auto kernelLabels = kernel->getExtraContributionLabels();
      if (static_cast<int>(kernelLabels.size()) !=
          kernel->getNbExtraContributions()) {
        throw std::runtime_error(fmt::format(
            "The number of labels returned by {} contains {} string,"
            "but it should contain {} (the number of contributions)",
            kernel->type(), kernelLabels.size(),
            kernel->getNbExtraContributions()));
      }
      labels.insert(labels.end(), kernelLabels.begin(), kernelLabels.end());
    }
    if (static_cast<int>(labels.size()) != m_nbContributions) {
      throw std::runtime_error(
          "The number of label is different from the number of contributions");
    }
    return labels;
  }
};
}  // namespace pairkernels

#endif
