// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#ifndef PBC_H
#define PBC_H

#include <cmath>

#include "util/array.h"
#include "util/moves.h"
#include "util/quaternions/quat.h"

namespace util {
namespace pbc {

/**
 * @brief ToClosest Return the value inPosition rounded between
 * [-inBoxSize/2;inBoxSize/2)
 * This method might be used to have the closest distance between two
 * beads regarding their possible images
 * @param inPosition a value possible greater/lower than (-)inBoxSize
 * @param inBoxSize the width of the box
 * @return a value in [-inBoxSize/2;inBoxSize/2)
 */
template <class RealType,
          typename = std::enable_if<std::is_floating_point<RealType>::value>>
inline RealType ToClosest(RealType inPosition, const RealType inBoxSize) {
  DEBUG_ASSERT(std::isfinite(inPosition),
               "ToClosest: inPosition is not finite");
  while (inPosition > inBoxSize * .5) {
    inPosition -= inBoxSize;
  }
  while (inPosition <= -inBoxSize * .5) {
    inPosition += inBoxSize;
  }
  return inPosition;
}

/**
 * @brief ToClosest apply ToClosest to the input vector
 * @param inPosition The distance between two points
 * @param inBoxSize The size of the box
 * @return The resulting vector where each value is in
 * [-inBoxSize/2;inBoxSize/2)
 */
template <class RealType,
          typename = std::enable_if<std::is_floating_point<RealType>::value>>
vec<RealType> ToClosest(const vec<RealType>& inPosition,
                        const vec<RealType>& inBoxSize) {
  return vec<RealType>(ToClosest<RealType>(inPosition[0], inBoxSize[0]),
                       ToClosest<RealType>(inPosition[1], inBoxSize[1]),
                       ToClosest<RealType>(inPosition[2], inBoxSize[2]));
}

/**
 * @brief ToClosest apply ToClosest to the input vector outPosition
 * @param outPosition The distance between two points
 * The resulting vector where each value is in [-inBoxSize/2;inBoxSize/2)
 * @param inBoxSize The size of the box
 */
template <class RealType,
          typename = std::enable_if<std::is_floating_point<RealType>::value>>
inline void ToClosest(vec<RealType>* const outPosition,
                      const vec<RealType>& inBoxSize) {
  (*outPosition)[0] = ToClosest<RealType>((*outPosition)[0], inBoxSize[0]);
  (*outPosition)[1] = ToClosest<RealType>((*outPosition)[1], inBoxSize[1]);
  (*outPosition)[2] = ToClosest<RealType>((*outPosition)[2], inBoxSize[2]);
}

/**
 * @brief ToBox Return the value inPosition rounded between [0;inBoxSize)
 * This method might be used to have map an image of a position in a grid
 * @param inPosition a value possible greater/lower than (-)inBoxSize
 * @param inBoxSize the width of the box
 * @return a value in [0;inBoxSize)
 */
template <class RealType,
          typename = std::enable_if<std::is_floating_point<RealType>::value>>
inline RealType ToBox(RealType position, const RealType boxWidth) {
  while (position > boxWidth) {
    position -= boxWidth;
  }
  while (position <= 0) {
    position += boxWidth;
  }
  return position;
}

/**
 * @brief ToBox apply ToBox to the input vector and return the result
 * @param inPosition a position in space
 * @param inBoxSize The size of the box
 * @return The resulting vector where each value is in [0;inBoxSize)
 */
template <class RealType,
          typename = std::enable_if<std::is_floating_point<RealType>::value>>
vec<RealType> ToBox(const vec<RealType>& inPosition,
                    const vec<RealType>& inBoxSize) {
  return util::vec<RealType>(ToBox<RealType>(inPosition[0], inBoxSize[0]),
                             ToBox<RealType>(inPosition[1], inBoxSize[1]),
                             ToBox<RealType>(inPosition[2], inBoxSize[2]));
}

/**
 * @brief ToBox apply ToBox to the input vector and update it
 * @param outPosition a position in space update to [0;inBoxSize)
 * @param inBoxSize The size of the box
 */
template <class RealType,
          typename = std::enable_if<std::is_floating_point<RealType>::value>>
inline void ToCBox(vec<RealType>* const outPosition,
                   const vec<RealType>& inBoxSize) {
  (*outPosition)[0] = ToBox<double>((*outPosition)[0], inBoxSize[0]);
  (*outPosition)[1] = ToBox<double>((*outPosition)[1], inBoxSize[1]);
  (*outPosition)[2] = ToBox<double>((*outPosition)[2], inBoxSize[2]);
}

/**
 * @brief DistSquare calculate the squared distance with PBC
 * @param inDistanceVec the distance between to point (a-b)
 * @param inBoxSize the size of the box
 * @return ToClosest(inDistanceVec).^2
 */
template <class RealType,
          typename = std::enable_if<std::is_floating_point<RealType>::value>>
inline RealType DistSquare(const util::vec<RealType>& inDistanceVec,
                           const util::vec<RealType>& inBoxSize) {
  const auto closestDistance = ToClosest(inDistanceVec, inBoxSize);

  const RealType distancePow2 =
      closestDistance[0] * closestDistance[0]     // x^2
      + closestDistance[1] * closestDistance[1]   // y^2
      + closestDistance[2] * closestDistance[2];  // z^2
  return distancePow2;
}

/**
 * @brief DistSquare calculate the squared distance with PBC
 * @param inPosition1 first position
 * @param inPosition2 second position
 * @param inSimulationBoxSize the size of the box
 * @param outShiftDom2 the shift applied to inPosition2 to be close to
 * inPosition1
 * @return ToClosest(inDistanceVec).^2
 */
inline double DistSquareBetweenPoints(const util::rvec& inPosition1,
                                      const util::rvec& inPosition2,
                                      const util::rvec inSimulationBoxSize,
                                      util::rvec* outShiftDom2) {
  util::rvec shiftPosition2 = inPosition2;
  for (int idx = 0; idx < 3; ++idx) {
    while (shiftPosition2[idx] - inPosition1[idx] >
           inSimulationBoxSize[idx] * .5) {
      shiftPosition2[idx] -= inSimulationBoxSize[idx];
      (*outShiftDom2)[idx] -= inSimulationBoxSize[idx];
    }
    while (shiftPosition2[idx] - inPosition1[idx] <=
           -inSimulationBoxSize[idx] * .5) {
      shiftPosition2[idx] += inSimulationBoxSize[idx];
      (*outShiftDom2)[idx] += inSimulationBoxSize[idx];
    }
  }
  return (inPosition1[0] - shiftPosition2[0]) *
             (inPosition1[0] - shiftPosition2[0]) +
         (inPosition1[1] - shiftPosition2[1]) *
             (inPosition1[1] - shiftPosition2[1]) +
         (inPosition1[2] - shiftPosition2[2]) *
             (inPosition1[2] - shiftPosition2[2]);
}

/**
 * @brief applyPBC translate an array of positions
 * in order to have the center of mass inside the box.
 * The center of mass must be in a neighboring image.
 * The resulting positions might be outside of the box.
 * @param xyz the positions to translate
 * @param box the simulation box
 * @return the translated array
 */
template <class RealType,
          typename = std::enable_if<std::is_floating_point<RealType>::value>>
void applyPBCInPlace(const vec<RealType>& box, Array<RealType>* xyz) {
  auto com = util::centroid(*xyz);
  auto diff = vec<RealType>();
  // This assumes that we are always in a neighboring image.
  for (auto j = 0; j < 3; ++j) {
    if (com[j] < 0) {
      diff[j] = box[j];
    } else if (com[j] > box[j]) {
      diff[j] = -box[j];
    } else {
      diff[j] = 0;
    }
  }
  for (auto i = 0; i < (*xyz).rows(); ++i) {
    for (auto j = 0; j < 3; ++j) {
      (*xyz)(i, j) = (*xyz)(i, j) + diff[j];
    }
  }
}

template <class RealType,
          typename = std::enable_if<std::is_floating_point<RealType>::value>>
Array<RealType> applyPBC(const Array<RealType>& xyz, const vec<RealType>& box) {
  auto res = xyz;
  applyPBCInPlace(box, &res);
  return res;
}

/**
 * @brief find closes image coordinates. Given a vector `a` and `b` find the
 * vector `b'` that is closest to a in all periodic boxes.
 *
 */
template <typename RealType>
vec<RealType> closestImage(const vec<RealType>& a, const vec<RealType>& b,
                           const vec<RealType>& box) {
  const auto box2 = vec<RealType>(0.5 * box[0], .5 * box[1], .5 * box[2]);
  auto res = b;
  for (auto i = 0; i < 3; ++i) {
    const auto diff = a[i] - b[i];
    if (diff < -box2[i]) {
      res[i] -= box[i];
    } else if (diff > box2[i]) {
      res[i] += box[i];
    }
  }
  return res;
}

}  // namespace pbc
}  // namespace util

#endif  // MOVES_H
