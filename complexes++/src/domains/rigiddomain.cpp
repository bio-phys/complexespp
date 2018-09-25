// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#include <array>
#include <cmath>
#include <fmt/format.h>
#include <random>

#include "constants.h"
#include "domains/rigiddomain.h"
#include "energy/energy.h"
#include "util/array.h"
#include "util/linalg.h"
#include "util/moves.h"
#include "util/pbc.h"
#include "util/util.h"

namespace q = util::quaternions;

namespace domains {

Rigid::Rigid(const std::string& inTypename, int typeId_, int id_,
             std::vector<Bead> beads_, std::vector<double> charges_,
             std::vector<BeadChainID> beadChainIDs_,
             const Connections& connections_, util::rvec trans_, double phi)
    : AbstractDomain(inTypename, typeId_, id_, beads_, charges_, beadChainIDs_,
                     connections_),
      m_phi(phi),
      m_trans(trans_) {}

MovedDomain Rigid::move(const util::rvec& box, util::RNGEngine& rng) const {
  auto dist = std::uniform_real_distribution<>{0, 1};
  auto tmp = *this;
  auto xyz_ = tmp.xyz();
  auto move_type = std::string();

  if (dist(rng) < 0.5) {
    const auto dx = util::randomVec(m_trans, rng);
    xyz_ = util::translateByConstant(xyz_, dx);
    move_type = "trans";
  } else {
    const auto rot = util::randomRotation<double>(m_phi, rng);
    const auto com = util::centroid(xyz_);
    xyz_ = util::rotateWithMat(xyz_, com, rot);
    move_type = "rot";
  }

  tmp.setXyz(util::pbc::applyPBC(xyz_, box));
  // do not do the conversion to a unique pointer before. Otherwise the type
  // system will srew us.
  return MovedDomain::Success(std::move(tmp), move_type);
}

std::unique_ptr<AbstractDomain> Rigid::copy() const {
  auto tmp = *this;
  return std::make_unique<Rigid>(std::move(tmp));
}
}  // namespace domains
