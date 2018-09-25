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
#ifndef UTIL_H
#define UTIL_H

#include <algorithm>
#include <cstdlib>
#include <iterator>
#include <vector>

#include "complexesconfig.h"

// Sometimes I don't want to use parameters provided by the overloaded API's.
// This will silence compiler warnings and result in a noop. Be careful if x is
// volatile.
#define UNUSED(x) (void)x;

////////////////////////////////////////////////////

#ifdef NDEBUG
#define DEBUG_ASSERT(condition, ...) ((void)0)
#else

#include <iostream>

#define DEBUG_ASSERT(condition, ...)                                           \
  if (!(condition)) {                                                          \
    fmt::print(std::cerr, "An assert has failed : {} \n", #condition);         \
    fmt::print(std::cerr, "\t In file : {}\n", __FILE__);                      \
    fmt::print(std::cerr, "\t At line : {}\n", __LINE__);                      \
    fmt::print(std::cerr, "\t Log : ");                                        \
    fmt::print(std::cerr, __VA_ARGS__);                                        \
    fmt::print(std::cerr, "\n");                                               \
    throw std::runtime_error("Bad Assert Exit");                               \
  }
#endif

////////////////////////////////////////////////////

#ifdef USE_TIMINGOUTPUT

/**
 * \def TIMEZONE
 * If enabled in the compilation configure,
 * TIMEZONE let time scope by wrapping ScopeEvent class.
 *
 * @code  {
 * @code     TIMEZONE("I want to measure from here");
 * @code     {
 * @code         TIMEZONE("I also want to measure this section");
 * @code
 * @code         for(int idx = 0 ; idx < 10 ; ++idx){
 * @code             TIMEZONE("Low level function A");
 * @code         }
 * @code     }{
 * @code         {
 * @code             TIMEZONE("This is another section");
 * @code         }
 * @code         TIMEZONE("Start another work here");
 * @code
 * @code         TIMEZONE("Low level function A");
 * @code     }
 * @code }
 * @code {
 * @code     TIMEZONE("I want to measure from here");
 * @code }
 * @code {
 * @code     TIMEZONE_MULTI_REF("A multi ref section");
 * @code }
 * @code {
 * @code     TIMEZONE_MULTI_REF("A multi ref section");
 * @code }
 */

#include "scopeevent.h"

extern util::EventManager GlobalEventManager;

#define TIMEZONE_Core_Merge(x, y) x##y
#define TIMEZONE_Core_Pre_Merge(x, y) TIMEZONE_Core_Merge(x, y)

#define TIMEZONE(NAME)                                                         \
  util::ScopeEvent TIMEZONE_Core_Pre_Merge(____TIMEZONE_AUTO_ID, __LINE__)(    \
      NAME, GlobalEventManager, ScopeEventUniqueKey);
#define TIMEZONE_MULTI_REF(NAME)                                               \
  util::ScopeEvent TIMEZONE_Core_Pre_Merge(____TIMEZONE_AUTO_ID, __LINE__)(    \
      NAME, GlobalEventManager, ScopeEventMultiRefKey);

#define TIMEZONE_OMP_INIT_PRETASK(VARNAME)                                     \
  auto VARNAME##core = GlobalEventManager.getCurrentThreadEvent();             \
  auto VARNAME = &VARNAME##core;
#define TIMEZONE_OMP_TASK(NAME, VARNAME)                                       \
  util::ScopeEvent TIMEZONE_Core_Pre_Merge(____TIMEZONE_AUTO_ID, __LINE__)(    \
      NAME, GlobalEventManager, ScopeEventUniqueKey, *VARNAME);
#define TIMEZONE_OMP_PRAGMA_TASK_KEY(VARNAME)                                  \
  shared(GlobalEventManager) firstprivate(VARNAME)

#define TIMEZONE_OMP_INIT_PREPARALLEL(NBTHREADS)                               \
  GlobalEventManager.startParallelRegion(NBTHREADS);

#else

#define TIMEZONE(NAME)
#define TIMEZONE_MULTI_REF(NAME)
#define TIMEZONE_OMP_INIT_PRETASK(VARNAME)
#define TIMEZONE_OMP_TASK(NAME, VARNAME)
#define TIMEZONE_OMP_PRAGMA_TASK_KEY(VARNAME)
#define TIMEZONE_OMP_INIT_PREPARALLEL(NBTHREADS)

#endif

////////////////////////////////////////////////////

namespace util {

template <typename Condition, typename T = void>
using EnableIf_t = typename std::enable_if<Condition::value, T>::type;

template <typename Iterator, typename IteratorTag>
using isSameIteratorCond =
    std::is_same<IteratorTag,
                 typename std::iterator_traits<Iterator>::iterator_category>;

template <typename Iter, typename T,
          typename = EnableIf_t<
              isSameIteratorCond<Iter, std::random_access_iterator_tag>>>
int indexOf(const Iter &begin, const Iter &end, const T &val) {
  const auto it = std::find(begin, end, val);
  if (end == it) {
    return -1;
  }
  return std::distance(begin, it);
}

template <class T> struct is_std_vector { static const int value = false; };

template <class T> struct is_std_vector<std::vector<T>> {
  static const int value = true;
};

} // namespace util
#endif // UTIL_H
