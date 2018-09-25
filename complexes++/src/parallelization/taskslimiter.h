// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#ifndef TASKSLIMITER_H
#define TASKSLIMITER_H

namespace vectorization {

class TasksLimiter {
  bool m_createTask;

 public:
  TasksLimiter() : m_createTask(true) {}

  void setEnableTasks(const bool inEnable) {
#pragma omp atomic write
    m_createTask = inEnable;
  }

  bool shouldCreateTasks() const {
    bool res;
#pragma omp atomic read
    res = m_createTask;
    return res;
  }

  static TasksLimiter Controller;
};
}

#endif
