// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#ifndef INTERACTIONBUILDER_HPP
#define INTERACTIONBUILDER_HPP

#include <memory>

#include "abstractinteractionalgorithm.h"
#include "cutoffgrid/codensegridcontainer.h"
#include "cutoffgrid/cosparsegridcontainer.h"
#include "cutoffgrid/cutoffinteractions.h"
#include "fullinteractions.h"

#include "domains/abstractdomain.h"
#include "io/io.h"
#include "setup/config.h"
#include "util/log.h"
#include "util/util.h"

/**
 * This function builds an interaction algorithm class
 * from the YAML config.
 */
template <class RealType>
std::unique_ptr<AbstractInteractionAlgorithm<RealType>>
buildInteractionComputer(std::shared_ptr<domains::Domains> doms,
                         const util::rvec& box, const setup::Config& conf) {
  TIMEZONE("buildInteractionComputer");
  const bool enableCutoff =
      conf.value<bool>("montecarlo.short-range-cutoff.enable");

  if (enableCutoff) {
    const RealType radius =
        conf.value<RealType>("montecarlo.short-range-cutoff.radius");
    const std::string cellContainers =
        conf.value<std::string>("montecarlo.short-range-cutoff.container");

    if (cellContainers == "dense") {
      util::Log("Config : cutoff - dense - radius = {} \n", radius);
      return std::make_unique<cutoffgrid::CutoffInteractions<
          RealType, cutoffgrid::CoDenseGridContainer>>(radius, box, doms);
    } else if (cellContainers == "sparse") {
      util::Log("Config : cutoff - sparse - radius = {} \n", radius);
      return std::make_unique<cutoffgrid::CutoffInteractions<
          RealType, cutoffgrid::CoSparseGridContainer>>(radius, box, doms);
    } else {
      throw std::invalid_argument(
          "key montecarlo.short-range-cutoff.container must either dense or "
          "sparse");
    }
  } else {
    util::Log("Config : full \n");
    return std::make_unique<FullInteractions<RealType>>(box, doms);
  }
}

#endif
