// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#ifndef FULLINTERACTIONS_HPP
#define FULLINTERACTIONS_HPP

#include <vector>

#include "parallelization/ompmanager.h"

#include "abstractinteractionalgorithm.h"
#include "domains/abstractdomain.h"
#include "energy/energymatrix.h"
#include "pairkernels/pairkernelmanager.h"
#include "util/log.h"
#include "util/pbc.h"

/**
 * This class provide the 2 interaction computation
 * methods. It does the computation in a full/direct way.
 * It just iterate over the domains and compute all the
 * interactions.
 * Therefor its complexity is high.
 */
template <class RealType>
class FullInteractions : public AbstractInteractionAlgorithm<RealType> {
 protected:
  static const int m_NbNeighbors = 26;

  const util::rvec m_boxSize;

  const std::shared_ptr<domains::Domains> m_allDomains;

 public:
  void serialize(io::Serializer& serializer) const final {
    AbstractInteractionAlgorithm<RealType>::serializeCore(serializer);
    serializer.append(m_boxSize, "m_boxSize");
  }

  FullInteractions(io::Deserializer& deserializer,
                   const std::shared_ptr<domains::Domains> inAllDomains)
      : AbstractInteractionAlgorithm<RealType>(deserializer),
        m_boxSize(deserializer.restore<decltype(m_boxSize)>("m_boxSize")),
        m_allDomains(inAllDomains) {}

  explicit FullInteractions(
      const util::rvec& inBoxSize,
      const std::shared_ptr<domains::Domains> inAllDomains)
      : m_boxSize(inBoxSize), m_allDomains(inAllDomains) {}

  void updateDomain(const int inDomainId) final { UNUSED(inDomainId); }

  void resetAllDomains(const util::rvec& newBox) final { UNUSED(newBox); }

  void computeAll(const util::rvec& box, const energy::ForceField& forcefield,
                  const pairkernels::PairKernelManager& inKernels,
                  energy::EnergyMatrix<RealType>& outRes) const final {
    TIMEZONE_OMP_INIT_PRETASK(timeZoneTaskKey);
    const int idxThreadInserted = omp_get_thread_num();

    bool shouldCreateTask =
        vectorization::TasksLimiter::Controller.shouldCreateTasks();
    omp_lock_t energyMatrixMutex;
    omp_init_lock(&energyMatrixMutex);

    std::vector<energy::EnergyMatrix<RealType>> connectionsBufferForThreads(
        omp_get_num_threads(), energy::EnergyMatrix<RealType>(0, 0, 0));

#pragma omp taskgroup
    {
      // The method iterates on all the domains in two loops
      // to compute the interactions between each pair.

      for (int idxDomTarget = 0;
           idxDomTarget < static_cast<int>(m_allDomains->size());
           ++idxDomTarget) {
        for (int idxDomSource = 0;
             idxDomSource < static_cast<int>(m_allDomains->size());
             ++idxDomSource) {
          if (shouldCreateTask == false) {
            shouldCreateTask =
                vectorization::TasksLimiter::Controller.shouldCreateTasks();
          }

          if (idxDomTarget != idxDomSource) {
#pragma omp task default(shared) firstprivate( \
    idxDomTarget, idxDomSource, idxThreadInserted) if (shouldCreateTask)
            {
              energy::EnergyMatrixBuffer<> bufferOutRes(outRes,
                                                        energyMatrixMutex);

              util::GlobalLog::redirectLog(idxThreadInserted);
              // Compute between idxDomTarget and idxDomSource
              inKernels
                  .getKernel((*m_allDomains)[idxDomTarget]->typeId(),
                             (*m_allDomains)[idxDomSource]->typeId())
                  .compute(bufferOutRes.getAccesser(
                               (*m_allDomains)[idxDomTarget]->id(),
                               (*m_allDomains)[idxDomSource]->id()),
                           *(*m_allDomains)[idxDomTarget],
                           *(*m_allDomains)[idxDomSource], box, forcefield);
            }
          }
        }
      }

      for (int idxDomTarget = 0;
           idxDomTarget < static_cast<int>(m_allDomains->size());
           ++idxDomTarget) {
// Add the connection energy of the idxDomTarget domain
#pragma omp task default(shared) firstprivate(idxDomTarget, idxThreadInserted) \
    TIMEZONE_OMP_PRAGMA_TASK_KEY(timeZoneTaskKey) if (shouldCreateTask)
        {
          util::GlobalLog::redirectLog(idxThreadInserted);
          TIMEZONE_OMP_TASK("connections", timeZoneTaskKey);
          if (connectionsBufferForThreads[omp_get_thread_num()]
                  .nbRowDomains() == 0) {
            connectionsBufferForThreads[omp_get_thread_num()] =
                energy::EnergyMatrix<RealType>(1UL, m_allDomains->size(), 1);
          }
          energy::EnergyMatrix<RealType>& connectionsBuffer =
              (connectionsBufferForThreads[omp_get_thread_num()]);
          connectionsBuffer.reset();
          (*m_allDomains)[idxDomTarget]->energyForAllConnections(
              (*m_allDomains), box, forcefield, connectionsBuffer);
          omp_set_lock(&energyMatrixMutex);
          for (int idxDomCon = 0; idxDomCon < connectionsBuffer.nbColDomains();
               ++idxDomCon) {
            outRes.addEnergyConnections(
                (*m_allDomains)[idxDomTarget]->id(), idxDomCon,
                connectionsBuffer.getEnergy(0, idxDomCon));
          }
          omp_unset_lock(&energyMatrixMutex);
        }
      }
    }
#if _OPENMP < 201307
    if (shouldCreateTask) {
#pragma omp taskwait
    }
#endif
  }

  void computeForOneDomain(const int inDomainId, const util::rvec& box,
                           const energy::ForceField& forcefield,
                           const pairkernels::PairKernelManager& inKernels,
                           energy::EnergyMatrix<RealType>& outRes) const final {
    DEBUG_ASSERT(inDomainId == (*m_allDomains)[inDomainId]->id(),
                 "The domain at position inDomainId has not an id equal to "
                 "inDomainId");
    DEBUG_ASSERT(outRes.nbRowDomains() == 1 &&
                     outRes.nbColDomains() == m_allDomains->size() &&
                     outRes.nbContributions() == inKernels.getNbContributions(),
                 "Invalid outRes size");
    TIMEZONE_OMP_INIT_PRETASK(timeZoneTaskKey);
    const int idxThreadInserted = omp_get_thread_num();
    omp_lock_t energyMatrixMutex;
    omp_init_lock(&energyMatrixMutex);

    // Here it is similar to the idxDomTarget but we do it only for the
    // inDomainId domain

    bool shouldCreateTask =
        vectorization::TasksLimiter::Controller.shouldCreateTasks();

    std::vector<std::unique_ptr<energy::EnergyMatrixBuffer<>>>
        energyBufferForThreads(omp_get_num_threads());
#pragma omp taskgroup
    {
      // Compute the energy with all other domains (including itself)
      for (int idxDomSource = 0;
           idxDomSource < static_cast<int>(m_allDomains->size());
           ++idxDomSource) {
        if (shouldCreateTask == false) {
          shouldCreateTask =
              vectorization::TasksLimiter::Controller.shouldCreateTasks();
        }

        if (inDomainId != idxDomSource) {
#pragma omp task default(shared) \
    firstprivate(idxDomSource, idxThreadInserted) if (shouldCreateTask)
          {
            util::GlobalLog::redirectLog(idxThreadInserted);
            energy::EnergyMatrixBuffer<> bufferOutRes(outRes,
                                                      energyMatrixMutex);
            inKernels
                .getKernel((*m_allDomains)[inDomainId]->typeId(),
                           (*m_allDomains)[idxDomSource]->typeId())
                .compute(bufferOutRes.getAccesser(
                             0, (*m_allDomains)[idxDomSource]->id()),
                         *(*m_allDomains)[inDomainId],
                         *(*m_allDomains)[idxDomSource], box, forcefield);
          }
        }
      }

      {
        omp_set_lock(&energyMatrixMutex);
        // Add the connections
        TIMEZONE_OMP_TASK("connections", timeZoneTaskKey);
        (*m_allDomains)[inDomainId]->energyForAllConnections(
            (*m_allDomains), box, forcefield, outRes);
        omp_unset_lock(&energyMatrixMutex);
      }
    }
#if _OPENMP < 201307
    if (shouldCreateTask) {
#pragma omp taskwait
    }
#endif
  }

  std::string type() const final { return FullInteractions::Type(); }

  static std::string Type() { return __PRETTY_FUNCTION__; }
};
REBUILDER_REGISTER(FullInteractions<double>);

#endif
