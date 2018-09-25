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
#ifndef CUTOFFINTERACTIONS_HPP
#define CUTOFFINTERACTIONS_HPP

#include "parallelization/ompmanager.h"

#include "cutoffgrid/cogrid.h"
#include "domains/abstractdomain.h"
#include "energy/energymatrix.h"
#include "interactions/abstractinteractionalgorithm.h"
#include "pairkernels/pairkernelmanager.h"
#include "parallelization/ompmanager.h"
#include "parallelization/taskslimiter.h"
#include "util/pbc.h"

namespace cutoffgrid {
/**
 * This class provide the 2 interaction computation
 * methods. It does the computation using a cutoff approach.
 * It has a low complexity but it has extra cost to keep
 * the m_grid up-to-date (in the CoGrid).
 */
template <class RealType, class GridContainerClass>
class CutoffInteractions : public AbstractInteractionAlgorithm<RealType> {
protected:
  static const int m_NbNeighbors = 26;

  const RealType m_coRadius;
  const util::rvec m_boxSize;
  //< The cutoff grid
  CoGrid<RealType, GridContainerClass> m_grid;

  const std::shared_ptr<domains::Domains> m_allDomains;

  size_t m_minMaxAvgNbCells[3];
  int m_countNbLoops;
  int m_nbTasksForOne;

public:
  virtual void serialize(io::Serializer &serializer) const final {
    AbstractInteractionAlgorithm<RealType>::serializeCore(serializer);
    serializer.append(m_coRadius, "m_coRadius");
    serializer.append(m_boxSize, "m_boxSize");
  }

  CutoffInteractions(io::Deserializer &deserializer,
                     const std::shared_ptr<domains::Domains> &inAllDomains)
      : AbstractInteractionAlgorithm<RealType>(deserializer),
        m_coRadius(deserializer.restore<decltype(m_coRadius)>("m_coRadius")),
        m_boxSize(deserializer.restore<decltype(m_boxSize)>("m_boxSize")),
        m_grid(m_coRadius, m_boxSize, inAllDomains->size()),
        m_allDomains(inAllDomains), m_nbTasksForOne(omp_get_max_threads()) {
    for (int idxDom = 0; idxDom < static_cast<int>(inAllDomains->size());
         ++idxDom) {
      m_grid.addDomain((*inAllDomains)[idxDom]);
    }
    m_minMaxAvgNbCells[0] = m_grid.getNbExistingCells();
    m_minMaxAvgNbCells[1] = m_grid.getNbExistingCells();
    m_minMaxAvgNbCells[2] = m_grid.getNbExistingCells();
    m_countNbLoops = 1;

    warnIfSlowCutoff();
  }

  explicit CutoffInteractions(
      const RealType inCoRadius, const util::rvec &inBoxSize,
      const std::shared_ptr<domains::Domains> inAllDomains)
      : m_coRadius(inCoRadius), m_boxSize(inBoxSize),
        m_grid(inCoRadius, inBoxSize, inAllDomains->size()),
        m_allDomains(inAllDomains), m_nbTasksForOne(omp_get_max_threads()) {
    for (int idxDom = 0; idxDom < static_cast<int>(inAllDomains->size());
         ++idxDom) {
      m_grid.addDomain((*inAllDomains)[idxDom]);
    }
    m_minMaxAvgNbCells[0] = m_grid.getNbExistingCells();
    m_minMaxAvgNbCells[1] = m_grid.getNbExistingCells();
    m_minMaxAvgNbCells[2] = m_grid.getNbExistingCells();
    m_countNbLoops = 1;

    warnIfSlowCutoff();
  }

  ~CutoffInteractions() {
    util::Log("@Min-nb-cells = {}s\n", m_minMaxAvgNbCells[0]);
    util::Log("@Max-nb-cells = {}s\n", m_minMaxAvgNbCells[1]);
    util::Log("@Avg-nb-cells = {}s\n", m_minMaxAvgNbCells[2]);
    util::Log("@Grid-size = {}s\n", m_grid.getNbMaximumCells());
  }

  void updateDomain(const int inDomainId) final {
    m_grid.updateDomain((*m_allDomains)[inDomainId]);
    m_minMaxAvgNbCells[0] =
        std::min(m_minMaxAvgNbCells[0], m_grid.getNbExistingCells());
    m_minMaxAvgNbCells[1] =
        std::max(m_minMaxAvgNbCells[1], m_grid.getNbExistingCells());
    m_minMaxAvgNbCells[2] = ((m_minMaxAvgNbCells[2] * m_countNbLoops) +
                             m_grid.getNbExistingCells()) /
                            (m_countNbLoops + 1);
    m_countNbLoops += 1;
  }

  void resetAllDomains(const util::rvec &newBox) final {
    m_grid = CoGrid<RealType, GridContainerClass>(m_coRadius, newBox,
                                                  m_allDomains->size());
    for (int idxDom = 0; idxDom < static_cast<int>(m_allDomains->size());
         ++idxDom) {
      m_grid.addDomain((*m_allDomains)[idxDom]);
    }
  }

  void computeAll(const util::rvec &box, const energy::ForceField &forcefield,
                  const pairkernels::PairKernelManager &inKernels,
                  energy::EnergyMatrix<RealType> &outRes) const final {
    std::vector<energy::EnergyMatrix<RealType>> connectionsBufferForThreads(
        omp_get_num_threads(), energy::EnergyMatrix<RealType>(0, 0, 0));

    const int idxThreadInserted = omp_get_thread_num();
    omp_lock_t energyMatrixMutex;
    omp_init_lock(&energyMatrixMutex);

    bool shouldCreateTask =
        vectorization::TasksLimiter::Controller.shouldCreateTasks();

    TIMEZONE_OMP_INIT_PRETASK(timeZoneTaskKey)
#pragma omp taskgroup
    {
      // Here we iterate over all the cells of the m_grid
      for (const CoCell &currentCellIter : m_grid) {
        const CoCell *currentCellPtr = &currentCellIter;

        if (shouldCreateTask == false) {
          shouldCreateTask =
              vectorization::TasksLimiter::Controller.shouldCreateTasks();
        }

#pragma omp task default(shared)                                               \
    firstprivate(currentCellPtr, idxThreadInserted) if (shouldCreateTask)
        {
          util::GlobalLog::redirectLog(idxThreadInserted);
          const CoCell &currentCell = *currentCellPtr;

          energy::EnergyMatrixBuffer<> bufferOutRes(outRes, energyMatrixMutex);

          // First we compute the interaction inside the cell with a double loop
          for (int idxDomTarget = 0; idxDomTarget < currentCell.getNbDomains();
               ++idxDomTarget) {
            const CoInterval &domainTarget =
                currentCell.getInterval(idxDomTarget);

            DEBUG_ASSERT(domainTarget.getDomainId() ==
                             (*m_allDomains)[domainTarget.getDomainId()]->id(),
                         "domainTarget.getDomainId() does not match");

            for (int idxDomSource = 0;
                 idxDomSource < currentCell.getNbDomains(); ++idxDomSource) {
              const CoInterval &domainSource =
                  currentCell.getInterval(idxDomSource);

              DEBUG_ASSERT(
                  domainSource.getDomainId() ==
                      (*m_allDomains)[domainSource.getDomainId()]->id(),
                  "domainSource.getDomainId() does not match");

              // Compute idxDomTarget idxDomSource interactions
              if (domainTarget.getDomainId() != domainSource.getDomainId()) {
                inKernels
                    .getKernel(
                        (*m_allDomains)[domainTarget.getDomainId()]->typeId(),
                        (*m_allDomains)[domainSource.getDomainId()]->typeId())
                    .compute(
                        bufferOutRes.getAccesser(domainTarget.getDomainId(),
                                                 domainSource.getDomainId()),
                        *(*m_allDomains)[domainTarget.getDomainId()],
                        *(*m_allDomains)[domainSource.getDomainId()], box,
                        forcefield, domainTarget.getInterval(),
                        domainSource.getInterval());
              }
            }
          }

          // Then between the current cell and the neigbhors
          std::array<const CoCell *, m_NbNeighbors> neighbors;
          const int nbNeighbors =
              m_grid.getPeriodicCellNeighbors(currentCell, &neighbors);

          for (int idxNeig = 0; idxNeig < nbNeighbors; ++idxNeig) {
            for (int idxDomTarget = 0;
                 idxDomTarget < currentCell.getNbDomains(); ++idxDomTarget) {
              const CoInterval &domainTarget =
                  currentCell.getInterval(idxDomTarget);

              DEBUG_ASSERT(
                  domainTarget.getDomainId() ==
                      (*m_allDomains)[domainTarget.getDomainId()]->id(),
                  "domainTarget.getDomainId() does not match");

              for (int idxDomSource = 0;
                   idxDomSource < neighbors[idxNeig]->getNbDomains();
                   ++idxDomSource) {
                const CoInterval &domainSource =
                    neighbors[idxNeig]->getInterval(idxDomSource);

                DEBUG_ASSERT(
                    domainSource.getDomainId() ==
                        (*m_allDomains)[domainSource.getDomainId()]->id(),
                    "domainSource.getDomainId() does not match");

                // Compute idxDomTarget idxDomSource (taken in a neighbor cell)
                // interactions
                if (domainTarget.getDomainId() != domainSource.getDomainId()) {
                  inKernels
                      .getKernel(
                          (*m_allDomains)[domainTarget.getDomainId()]->typeId(),
                          (*m_allDomains)[domainSource.getDomainId()]->typeId())
                      .compute(
                          bufferOutRes.getAccesser(domainTarget.getDomainId(),
                                                   domainSource.getDomainId()),
                          *(*m_allDomains)[domainTarget.getDomainId()],
                          *(*m_allDomains)[domainSource.getDomainId()], box,
                          forcefield, domainTarget.getInterval(),
                          domainSource.getInterval());
                }
              }
            }
          }
        }
      }
    }

#if _OPENMP < 201307
    if (shouldCreateTask) {
#pragma omp taskwait
    }
#endif

#pragma omp taskgroup
    {
      // Finally we add the connections and inner energy
      for (int idxDomTarget = 0;
           idxDomTarget < static_cast<int>(m_allDomains->size());
           ++idxDomTarget) {
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
          energy::EnergyMatrix<RealType> &connectionsBuffer =
              (connectionsBufferForThreads[omp_get_thread_num()]);
          connectionsBuffer.reset();
          (*m_allDomains)[idxDomTarget]->energyForAllConnections(
              (*m_allDomains), box, forcefield, connectionsBuffer);
          omp_set_lock(&energyMatrixMutex);
          for (int idxDomCon = 0;
               idxDomCon < static_cast<int>(connectionsBuffer.nbColDomains());
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

    omp_destroy_lock(&energyMatrixMutex);
  }

  void computeForOneDomain(const int inDomainId, const util::rvec &box,
                           const energy::ForceField &forcefield,
                           const pairkernels::PairKernelManager &inKernels,
                           energy::EnergyMatrix<RealType> &outRes) const final {
    DEBUG_ASSERT(inDomainId == (*m_allDomains)[inDomainId]->id(),
                 "The domain at position inDomainId has not an id equal to "
                 "inDomainId");
    DEBUG_ASSERT(outRes.nbRowDomains() == 1 &&
                     outRes.nbColDomains() == m_allDomains->size() &&
                     outRes.nbContributions() == inKernels.getNbContributions(),
                 "Invalid outRes size");
    // This function is similar to computeAll except that
    // there is only one target

    // We take the list of cell of the target domain
    const std::vector<CoDomainCellLink> &domainCells =
        m_grid.getDomainCells(inDomainId);

    bool shouldCreateTask =
        vectorization::TasksLimiter::Controller.shouldCreateTasks();
    omp_lock_t energyMatrixMutex;
    omp_init_lock(&energyMatrixMutex);

    const int idxThreadInserted = omp_get_thread_num();

    TIMEZONE_OMP_INIT_PRETASK(timeZoneTaskKey)
#pragma omp taskgroup
    {
      int nbCellsPerTask = 1;
      if (shouldCreateTask) {
        nbCellsPerTask =
            std::max(1, static_cast<int>(domainCells.size()) / m_nbTasksForOne);
      }

      for (int idxCellLink = 0;
           idxCellLink < static_cast<int>(domainCells.size());
           idxCellLink += nbCellsPerTask) {
        if (shouldCreateTask == false) {
          shouldCreateTask =
              vectorization::TasksLimiter::Controller.shouldCreateTasks();

          if (shouldCreateTask) {
            nbCellsPerTask = std::max(
                1, (static_cast<int>(domainCells.size()) - idxCellLink) /
                       m_nbTasksForOne);
          }
        }

#pragma omp task default(shared)                                               \
    firstprivate(idxCellLink, nbCellsPerTask, idxThreadInserted)               \
        TIMEZONE_OMP_PRAGMA_TASK_KEY(timeZoneTaskKey) if (shouldCreateTask)
        {
          util::GlobalLog::redirectLog(idxThreadInserted);
          TIMEZONE_OMP_TASK("compute", timeZoneTaskKey);
          const int endIdxCell = std::min(static_cast<int>(domainCells.size()),
                                          idxCellLink + nbCellsPerTask);
          energy::EnergyMatrixBuffer<> bufferOutRes(outRes, energyMatrixMutex);

          for (int idxCellTask = idxCellLink; idxCellTask < endIdxCell;
               ++idxCellTask) {
            const CoDomainCellLink &cellLink = domainCells[idxCellTask];
            const long long cellIndex = cellLink.getCellIndex();
            const CoCell *cellPtr = &m_grid.getCell(cellIndex);
            const CoInterval *domainTargetPtr =
                &cellPtr->getInterval(cellLink.getInsertPosInList());

            DEBUG_ASSERT(
                domainTargetPtr->getDomainId() == inDomainId,
                "The domain from domainTarget should be equal to inDomainId");
            {
              const CoCell &cell = *cellPtr;
              const CoInterval &domainTarget = *domainTargetPtr;
              // First the interaction inside the cell
              for (int idxDom = 0; idxDom < cell.getNbDomains(); ++idxDom) {
                const CoInterval &domainSource = cell.getInterval(idxDom);

                DEBUG_ASSERT(
                    domainSource.getDomainId() ==
                        (*m_allDomains)[domainSource.getDomainId()]->id(),
                    "domainSource.getDomainId() does not match");
                if (inDomainId != domainSource.getDomainId()) {
                  inKernels
                      .getKernel(
                          (*m_allDomains)[inDomainId]->typeId(),
                          (*m_allDomains)[domainSource.getDomainId()]->typeId())
                      .compute(bufferOutRes.getAccesser(
                                   0, domainSource.getDomainId()),
                               *(*m_allDomains)[inDomainId],
                               *(*m_allDomains)[domainSource.getDomainId()],
                               box, forcefield, domainTarget.getInterval(),
                               domainSource.getInterval());
                }
              }
            }

            {
              const CoCell &cell = *cellPtr;
              const CoInterval &domainTarget = *domainTargetPtr;
              // The the interactions with the neighbors
              std::array<const CoCell *, m_NbNeighbors> neighbors;
              const int nbNeighbors =
                  m_grid.getPeriodicCellNeighbors(cell, &neighbors);

              for (int idxNeig = 0; idxNeig < nbNeighbors; ++idxNeig) {
                for (int idxDom = 0;
                     idxDom < neighbors[idxNeig]->getNbDomains(); ++idxDom) {
                  const CoInterval &domainSource =
                      neighbors[idxNeig]->getInterval(idxDom);

                  DEBUG_ASSERT(
                      domainSource.getDomainId() ==
                          (*m_allDomains)[domainSource.getDomainId()]->id(),
                      "domainSource.getDomainId() does not match");

                  if (inDomainId != domainSource.getDomainId()) {
                    inKernels
                        .getKernel((*m_allDomains)[inDomainId]->typeId(),
                                   (*m_allDomains)[domainSource.getDomainId()]
                                       ->typeId())
                        .compute(bufferOutRes.getAccesser(
                                     0, domainSource.getDomainId()),
                                 *(*m_allDomains)[inDomainId],
                                 *(*m_allDomains)[domainSource.getDomainId()],
                                 box, forcefield, domainTarget.getInterval(),
                                 domainSource.getInterval());
                  }
                }
              }
            }
          }
        }
      }

      {
        // Add the connections
        TIMEZONE_OMP_TASK("connections", timeZoneTaskKey);
        omp_set_lock(&energyMatrixMutex);
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

    omp_destroy_lock(&energyMatrixMutex);
  }

  std::string type() const final { return CutoffInteractions::Type(); }

  static std::string Type() { return __PRETTY_FUNCTION__; }

  void warnIfSlowCutoff() const {
    // Count the total number of full interactions (approximatly!)
    size_t totalNbBeads = 0;
    for (int idxDom = 0; idxDom < static_cast<int>(m_allDomains->size());
         ++idxDom) {
      totalNbBeads += (*m_allDomains)[idxDom]->nBeads();
    }
    const size_t totalNbFullInteractions = totalNbBeads * totalNbBeads;

    // Count the total number of cutoff interactions
    size_t totalNbCoInteractions = 0;
    for (const CoCell &currentCell : m_grid) {
      // First we compute the interaction inside the cell with a double loop
      for (int idxDomTarget = 0; idxDomTarget < currentCell.getNbDomains();
           ++idxDomTarget) {
        const CoInterval &domainTarget = currentCell.getInterval(idxDomTarget);

        for (int idxDomSource = 0; idxDomSource < currentCell.getNbDomains();
             ++idxDomSource) {
          const CoInterval &domainSource =
              currentCell.getInterval(idxDomSource);

          totalNbCoInteractions += domainTarget.getNbElementsInInterval() *
                                   domainSource.getNbElementsInInterval();
        }
      }

      // Then between the current cell and the neigbhors
      std::array<const CoCell *, m_NbNeighbors> neighbors;
      const int nbNeighbors =
          m_grid.getPeriodicCellNeighbors(currentCell, &neighbors);

      for (int idxNeig = 0; idxNeig < nbNeighbors; ++idxNeig) {
        for (int idxDomTarget = 0; idxDomTarget < currentCell.getNbDomains();
             ++idxDomTarget) {
          const CoInterval &domainTarget =
              currentCell.getInterval(idxDomTarget);

          for (int idxDomSource = 0;
               idxDomSource < neighbors[idxNeig]->getNbDomains();
               ++idxDomSource) {
            const CoInterval &domainSource =
                neighbors[idxNeig]->getInterval(idxDomSource);

            totalNbCoInteractions += domainTarget.getNbElementsInInterval() *
                                     domainSource.getNbElementsInInterval();
          }
        }
      }
    }

    // Print warning if the number of interactions from cutoff is 80%
    // or higher than the number of interactions from the full method
    if (totalNbCoInteractions >= totalNbFullInteractions * 0.8) {
      util::Log("Warning -- From the given configuration, the Cutoff method\n");
      util::Log(
          "Warning -- might be slower than the direct/full computation.\n");
      util::Log(
          "Warning -- It is recommanded to test both before executing a long "
          "run.\n");
      util::Log(
          "Warning -- Cutoff interactions approximatly equal to {}% of the\n",
          100. * static_cast<double>(totalNbCoInteractions) /
              static_cast<double>(totalNbFullInteractions));
      util::Log(
          "Warning -- direct/full interactions, But the overhead of using a "
          "cutoff might result in longer runtime.\n");
    }
  }
};
using CutoffInteractions_double_cutoffgrid_CoDenseGridContainer =
    CutoffInteractions<double, cutoffgrid::CoDenseGridContainer>;
REBUILDER_REGISTER(CutoffInteractions_double_cutoffgrid_CoDenseGridContainer);
using CutoffInteractions_double_cutoffgrid_CoSparseGridContainer =
    CutoffInteractions<double, cutoffgrid::CoSparseGridContainer>;
REBUILDER_REGISTER(CutoffInteractions_double_cutoffgrid_CoSparseGridContainer);

} // namespace cutoffgrid
#endif // CUTOFFINTERACTIONS_HPP
