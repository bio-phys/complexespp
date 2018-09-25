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
} // namespace vectorization

#endif
