#ifndef ENERGY_PAIRKERNELBUILDER_H
#define ENERGY_PAIRKERNELBUILDER_H

#include <functional>
#include <sstream>
#include <unordered_map>

#include "pairkernels/ljhcutoffpairkernel.h"
#include "pairkernels/ljhpairkernel.h"
#include "pairkernels/membranepairkernel.h"
#include "pairkernels/nonepairkernel.h"
#include "pairkernels/repulsivepairkernel.h"
#include "pairkernels/softcore.h"
#include "pairkernels/wcapairkernel.h"

namespace pairkernels {

class PairKernelBuilder {
  using BuildFuncType = std::function<std::unique_ptr<AbstractPairKernel>()>;

  // To store extra kernels (type name => build function)
  static std::unordered_map<std::string, BuildFuncType> m_ExtraBuilders;

 public:
  // Build a kernel from its name,
  // The function first test the kernels provided by complexes,
  // before looking at the one registered with RegisterExtraBuild.
  static inline std::unique_ptr<AbstractPairKernel> Build(
      const std::string& kernelName) {
    if (kernelName == pairkernels::LJHPairKernel::GetName()) {
      return pairkernels::LJHPairKernel::BuildLJHPairKernel();
    } else if (kernelName == pairkernels::LJHCutoffPairKernel::GetName()) {
      return pairkernels::LJHCutoffPairKernel::BuildLJHCutoffPairKernel();
    } else if (kernelName == pairkernels::NonePairKernel::GetName()) {
      return pairkernels::NonePairKernel::BuildNonePairKernel();
    } else if (kernelName == pairkernels::RepulsivePairKernel::GetName()) {
      return pairkernels::RepulsivePairKernel::BuildRepulsivePairKernel();
    } else if (kernelName == pairkernels::MembranePairKernel::GetName()) {
      return pairkernels::MembranePairKernel::BuildMembranePairKernel();
    } else if (kernelName == pairkernels::WCAPairKernel::GetName()) {
      return pairkernels::WCAPairKernel::BuildWCAPairKernel();
    } else if (kernelName == pairkernels::SoftcorePairKernel::GetName()) {
      return pairkernels::SoftcorePairKernel::BuildSoftcorePairKernel();
    } else {
      std::unique_ptr<AbstractPairKernel> fromExtra =
          BuildFromExtra(kernelName);
      if (fromExtra) {
        return fromExtra;
      }

      std::stringstream strKernels;
      strKernels << pairkernels::LJHPairKernel::GetName() << ", "
                 << pairkernels::LJHCutoffPairKernel::GetName() << ", "
                 << pairkernels::NonePairKernel::GetName() << ", "
                 << pairkernels::WCAPairKernel::GetName() << ", "
                 << pairkernels::SoftcorePairKernel::GetName() << ", "
                 << pairkernels::RepulsivePairKernel::GetName();
      for (const auto& iter : m_ExtraBuilders) {
        strKernels << ", " << iter.first;
      }

      throw std::runtime_error(
          fmt::format("BuildPairKernel -- Pair kernel {} cannot be found, "
                      "available kernels are:\n{}",
                      kernelName, strKernels.str()));
    }
  }

  // Build a kernel from its name if it has been registered with
  // RegisterExtraBuild
  // It is useful for test kernels.
  static inline std::unique_ptr<AbstractPairKernel> BuildFromExtra(
      const std::string& kernelName) {
    for (const auto& iter : m_ExtraBuilders) {
      if (iter.first == kernelName) {
        return iter.second();
      }
    }

    return std::unique_ptr<AbstractPairKernel>();
  }

  // Register an extra kernel
  static inline bool RegisterExtraBuild(const std::string& kernelName,
                                        BuildFuncType buildFunc) {
    m_ExtraBuilders[kernelName] = buildFunc;
    return true;
  }
};
}  // namespace pairkernels

#endif
