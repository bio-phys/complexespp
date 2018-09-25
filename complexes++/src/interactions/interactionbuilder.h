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
                         const util::rvec &box, const setup::Config &conf) {
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
