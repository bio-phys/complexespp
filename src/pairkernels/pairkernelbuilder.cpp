#include "pairkernelbuilder.h"

std::unordered_map<std::string, pairkernels::PairKernelBuilder::BuildFuncType>
    pairkernels::PairKernelBuilder::m_ExtraBuilders;
