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
#ifndef TIMER_H
#define TIMER_H

#include <chrono>

namespace util {

/**
 * @file
 * The Timer class allow to measure elapsed time easily.
 * But to profile code and have the appropriate output
 * it is recommended to use ScopeEvent class.
 *
 * Each section to measure should be embraced by start/stop.
 * The measured time is given by "getElapsed".
 * The total time measured by a timer is given by "getCumulated".
 * Example :
 * @code Timer tm; // Implicit start
 * @code ...
 * @code tm.stop(); // stop the timer
 * @code tm.getElapsed(); // return the duration in s [A]
 * @code tm.start(); // restart the timer
 * @code ...
 * @code tm.stopAndGetElapsed(); // stop the timer and return the duraction in s
 * [B]
 * @code tm.getCumulated(); // Equal [A] + [B]
 */
class Timer {
  using double_second_time = std::chrono::duration<double, std::ratio<1, 1>>;

  std::chrono::high_resolution_clock::time_point
      m_start; ///< m_start time (start)
  std::chrono::high_resolution_clock::time_point m_end; ///< stop time (stop)
  std::chrono::nanoseconds m_cumulate;                  ///< the m_cumulate time

public:
  /// Constructor
  Timer() { start(); }

  /// Copy constructor
  Timer(const Timer &other) = delete;
  /// Copies an other timer
  Timer &operator=(const Timer &other) = delete;
  /// Move constructor
  Timer(Timer &&other) = default;
  /// Copies an other timer
  Timer &operator=(Timer &&other) = default;

  /** Rest all the values, and apply start */
  void reset() {
    m_start = std::chrono::high_resolution_clock::time_point();
    m_end = std::chrono::high_resolution_clock::time_point();
    m_cumulate = std::chrono::nanoseconds();
    start();
  }

  /** Start the timer */
  void start() { m_start = std::chrono::high_resolution_clock::now(); }

  /** Stop the current timer */
  void stop() {
    m_end = std::chrono::high_resolution_clock::now();
    m_cumulate +=
        std::chrono::duration_cast<std::chrono::nanoseconds>(m_end - m_start);
  }

  /** Return the elapsed time between start and stop (in second) */
  double getElapsed() const {
    return std::chrono::duration_cast<double_second_time>(
               std::chrono::duration_cast<std::chrono::nanoseconds>(m_end -
                                                                    m_start))
        .count();
  }

  /** Return the total counted time */
  double getCumulated() const {
    return std::chrono::duration_cast<double_second_time>(m_cumulate).count();
  }

  /** End the current counter (stop) and return the elapsed time */
  double stopAndGetElapsed() {
    stop();
    return getElapsed();
  }
};
} // namespace util

#endif
