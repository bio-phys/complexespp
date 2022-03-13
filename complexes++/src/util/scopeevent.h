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
#ifndef SCOPEEVENT_H
#define SCOPEEVENT_H
#define FMT_HEADER_ONLY
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <iostream>
#include <list>
#include <stack>
#include <string>
#include <unordered_map>
#include <vector>

#include "parallelization/ompmanager.h"
#include "timer.h"
#include "util.h"

namespace util {

//< To add it as friend of EventManager
class ScopeEvent;

/**
 * The EventManager class is in charge of the record of the
 * events and the hierarchy.
 * It is also responsible of the print out the records of the
 * managed events.
 * The ScopeEvent are linked by a EventManager, but
 * several EventManager can be used in an application.
 *
 * Example:
 * @code util::EventManager manager("Application name", std::cout);
 * @code
 * @code{
 * @code
 * @code    util::ScopeEvent event("I want to measure from here", manager,
 *ScopeEventUniqueKey);
 * @code    {
 * @code        util::ScopeEvent event1("I also want to measure this section",
 *manager, ScopeEventUniqueKey);
 * @code
 * @code        for(int idx = 0 ; idx < 10 ; ++idx){
 * @code            util::ScopeEvent eventA("Low level function A", manager,
 *ScopeEventUniqueKey);
 * @code        }
 * @code    }{
 * @code        {
 * @code            util::ScopeEvent event1("This is another section", manager,
 *ScopeEventUniqueKey);
 * @code        }
 * @code        util::ScopeEvent event2("Start another work here", manager,
 *ScopeEventUniqueKey);
 * @code
 * @code        util::ScopeEvent eventA("Low level function A", manager,
 *ScopeEventUniqueKey);
 * @code    }
 * @code}
 * @code{
 * @code    util::ScopeEvent eventBis("I want to measure from here", manager,
 *ScopeEventUniqueKey);
 * @code}
 *
 * It is recommended to have a look to USE_TIMINGOUTPUT/TIMEZONE to know more
 * about the default timing system.
 *
 * The default output is:
 * "@[KEY] = [TIME]s"
 *
 * If an event appear more than one (with the same stack) its output is:
 * "@[KEY] = [TIME]s (Min = [TIME]s ; Max = [TIME]s ; Average = [TIME] ;
 *Occurrence
 *= [TIME])"
 */
class EventManager {
protected:
  /**
   * The CoreEvent class represent an event.
   * It stores the duration/min/max/occurrence
   * and it has a stack representation of its parents
   * in order to distinct different call to the same event loger.
   */
  class CoreEvent {
  protected:
    //< Name of the event (from the user)
    const std::string m_name;
    //< Previous events (stack of parents)
    std::stack<CoreEvent *> m_parentStack;
    //< Current event children
    std::vector<CoreEvent *> m_children;

    //< Total execution time
    double m_totalTime;
    //< Minimum execution time
    double m_minTime;
    //< Maximum execution time
    double m_maxTime;
    //< Number of occurrence for this event
    int m_occurrence;
    //< Number of occurrence that are tasks for this event
    int m_nbTasks;
    //< Children lock
    omp_lock_t m_childrenLock;
    //< Children lock
    omp_lock_t m_updateLock;

  public:
    /** Create a core-event from the name and the current stack */
    CoreEvent(const std::string &inName,
              const std::stack<CoreEvent *> &inParentStack)
        : m_name(inName), m_parentStack(inParentStack), m_totalTime(0),
          m_minTime(std::numeric_limits<double>::max()),
          m_maxTime(std::numeric_limits<double>::min()), m_occurrence(0),
          m_nbTasks(0) {
      omp_init_lock(&m_childrenLock);
      omp_init_lock(&m_updateLock);
    }

    ~CoreEvent() {
      omp_destroy_lock(&m_childrenLock);
      omp_destroy_lock(&m_updateLock);
    }

    /** Add a record */
    void addRecord(const double inDuration, const bool isTask) {
#pragma omp atomic update
      m_totalTime += inDuration;
#pragma omp atomic update
      m_occurrence += 1;
#pragma omp flush // (m_minTime, m_maxTime)
      if (inDuration < m_minTime || m_maxTime < inDuration) {
        omp_set_lock(&m_updateLock);
        m_minTime = std::min(m_minTime, inDuration);
        m_maxTime = std::max(m_maxTime, inDuration);
        omp_unset_lock(&m_updateLock);
      }
      if (isTask) {
#pragma omp atomic update
        m_nbTasks += 1;
      }
    }

    const std::stack<CoreEvent *> &getParents() const { return m_parentStack; }

    std::stack<CoreEvent *> &getParents() { return m_parentStack; }

    void addChild(CoreEvent *inChild) {
      omp_set_lock(&m_childrenLock);
      m_children.push_back(inChild);
      omp_unset_lock(&m_childrenLock);
    }

    //! Must not be called during a paralle execution
    const std::vector<CoreEvent *> &getChildren() const {
      DEBUG_ASSERT(omp_in_parallel() == 0,
                   "getChildren cannot be called from a parrallel region");
      return m_children;
    }

    const std::string &getName() const { return m_name; }

    double getMin() const { return m_minTime; }

    double getMax() const { return m_maxTime; }

    int getOccurrence() const { return m_occurrence; }

    double getAverage() const {
      return m_totalTime / static_cast<double>(m_occurrence);
    }

    double getDuration() const { return m_totalTime; }

    int getNbTasks() const { return m_nbTasks; }
  };

  ///////////////////////////////////////////////////////////////

  //< The main node
  std::unique_ptr<CoreEvent> m_root;
  //< Output stream to print out
  std::ostream &m_outputStream;

  //< Current stack, there are one stack of stack per thread
  std::vector<std::stack<std::stack<CoreEvent *>>>
      m_currentEventsStackPerThread;
  //< All recorded events (that will then be delete at the end)
  std::unordered_multimap<std::string, CoreEvent *> m_records;
  //< Lock for m_records
  omp_lock_t m_recordsLock;

  /** Find a event from its name. If such even does not exist
   * the function creates one. If an event with the same name exists
   * but with a different stack, a new one is created.
   * It pushes the returned event in the stack.
   */
  CoreEvent *getEvent(const std::string &inName,
                      const std::string &inUniqueKey) {
    const std::string completeName = inName + inUniqueKey;
    CoreEvent *foundEvent = nullptr;

    omp_set_lock(&m_recordsLock);
    // find all events with this name
    auto range = m_records.equal_range(completeName);
    for (auto iter = range.first; iter != range.second; ++iter) {
      // events are equal if same name and same parents
      if ((*iter).second->getParents() ==
          m_currentEventsStackPerThread[omp_get_thread_num()].top()) {
        foundEvent = (*iter).second;
        break;
      }
    }

    // Keep the lock to ensure that not two threads create the same event

    if (!foundEvent) {
      // create this event
      foundEvent = new CoreEvent(
          inName, m_currentEventsStackPerThread[omp_get_thread_num()].top());
      m_currentEventsStackPerThread[omp_get_thread_num()].top().top()->addChild(
          foundEvent);
      m_records.insert({completeName, foundEvent});
    }
    omp_unset_lock(&m_recordsLock);

    m_currentEventsStackPerThread[omp_get_thread_num()].top().push(foundEvent);
    return foundEvent;
  }

  CoreEvent *getEventFromContext(const std::string &inName,
                                 const std::string &inUniqueKey,
                                 const std::stack<CoreEvent *> &inParentStack) {
    m_currentEventsStackPerThread[omp_get_thread_num()].push(inParentStack);
    return getEvent(inName, inUniqueKey);
  }

  /** Pop current event */
  void popEvent(const CoreEvent *eventToRemove) {
    DEBUG_ASSERT(
        m_currentEventsStackPerThread[omp_get_thread_num()].top().size() > 1,
        "Poped to many events, root event cannot be poped");
    // Comparing address is cheaper
    if (m_currentEventsStackPerThread[omp_get_thread_num()].top().top() !=
        eventToRemove) {
      throw std::runtime_error(
          "You must end events (ScopeEvent/TIMEZONE) in order.\n"
          "Please make sure that you only ask to the last event to finish.");
    }
    m_currentEventsStackPerThread[omp_get_thread_num()].top().pop();
  }

  /** Pop current context */
  void popContext(const CoreEvent *eventToRemove) {
    DEBUG_ASSERT(m_currentEventsStackPerThread[omp_get_thread_num()].size() > 1,
                 "Poped to many context");
    DEBUG_ASSERT(
        m_currentEventsStackPerThread[omp_get_thread_num()].top().size() > 1,
        "Poped to many events, root event cannot be poped");
    // Comparing address is cheaper
    if (m_currentEventsStackPerThread[omp_get_thread_num()].top().top() !=
        eventToRemove) {
      throw std::runtime_error(
          "You must end events (ScopeEvent/TIMEZONE) in order.\n"
          "Please make sure that you only ask to the last event to finish.");
    }
    m_currentEventsStackPerThread[omp_get_thread_num()].pop();
  }

public:
  /** Create an event manager */
  EventManager(const std::string &inAppName, std::ostream &inOutputStream)
      : m_root(new CoreEvent(inAppName, std::stack<CoreEvent *>())),
        m_outputStream(inOutputStream), m_currentEventsStackPerThread(1) {
    m_currentEventsStackPerThread[0].emplace();
    m_currentEventsStackPerThread[0].top().push(m_root.get());
    omp_init_lock(&m_recordsLock);
  }

  ~EventManager() throw() {
    show();

    DEBUG_ASSERT(m_currentEventsStackPerThread[0].size() == 1,
                 "Oups, the event-stack is corrupted, should 1 and is {}",
                 m_currentEventsStackPerThread[0].size());

    DEBUG_ASSERT(m_currentEventsStackPerThread[0].top().size() == 1,
                 "Oups, the event-stack is corrupted, should 1 and is {}",
                 m_currentEventsStackPerThread[0].top().size());

    omp_destroy_lock(&m_recordsLock);

    for (auto event : m_records) {
      delete event.second;
    }
  }

  void startParallelRegion(const int inNbThreads) {
    m_currentEventsStackPerThread.resize(1);
    m_currentEventsStackPerThread.resize(inNbThreads,
                                         m_currentEventsStackPerThread[0]);
  }

  void show() { show(m_outputStream); }

  void show(std::ostream &inOutputStream) const {
    std::stack<std::pair<int, const CoreEvent *>> events;

    for (int idx = static_cast<int>(m_root->getChildren().size()) - 1; idx >= 0;
         --idx) {
      events.push({0, m_root->getChildren()[idx]});
    }

    fmt::print(inOutputStream, "[TIMING] {}:\n", m_root->getName());

    while (events.size()) {
      const std::pair<int, const CoreEvent *> eventToShow = events.top();
      events.pop();

      int offsetTab = eventToShow.first;
      while (offsetTab--) {
        fmt::print(inOutputStream, "\t");
      }
      fmt::print(inOutputStream, "@{} = {}s", eventToShow.second->getName(),
                 eventToShow.second->getDuration());
      if (eventToShow.second->getOccurrence() != 1) {
        fmt::print(inOutputStream,
                   " (Min = {}s ; Max = {}s ; Average = {}s ; Occurrence = {})",
                   eventToShow.second->getMin(), eventToShow.second->getMax(),
                   eventToShow.second->getAverage(),
                   eventToShow.second->getOccurrence());
      }
      if (eventToShow.second->getNbTasks()) {
        fmt::print(inOutputStream, " => {} events are tasks",
                   eventToShow.second->getNbTasks());
      }

      fmt::print(inOutputStream, "\n");
      for (int idx =
               static_cast<int>(eventToShow.second->getChildren().size()) - 1;
           idx >= 0; --idx) {
        events.push(
            {eventToShow.first + 1, eventToShow.second->getChildren()[idx]});
      }
    }
  }

  std::stack<CoreEvent *> getCurrentThreadEvent() const {
    return m_currentEventsStackPerThread[omp_get_thread_num()].top();
  }

  friend ScopeEvent;
};

///////////////////////////////////////////////////////////////

/** A scope event should be used
 * to record the duration of a part of the code
 * (section, scope, etc.).
 * The timer is stoped automatically when the object is destroyed
 * or when "finish" is explicitely called.
 * The object cannot be copied/moved to ensure coherency in the
 * events hierarchy.
 */
class ScopeEvent {
protected:
  //< The manager to refer to
  EventManager &m_manager;
  //< The core event
  EventManager::CoreEvent *const m_event;
  //< Time to get elapsed time
  Timer m_timer;
  //< Is true if it has been created for task
  bool m_isTask;

public:
  ScopeEvent(const std::string &inName, EventManager &inManager,
             const std::string &inUniqueKey)
      : m_manager(inManager), m_event(inManager.getEvent(inName, inUniqueKey)),
        m_isTask(false) {
    m_timer.start();
  }

  ScopeEvent(const std::string &inName, EventManager &inManager,
             const std::string &inUniqueKey,
             const std::stack<EventManager::CoreEvent *> &inParentStack)
      : m_manager(inManager), m_event(inManager.getEventFromContext(
                                  inName, inUniqueKey, inParentStack)),
        m_isTask(true) {
    m_timer.start();
  }

  ~ScopeEvent() {
    m_event->addRecord(m_timer.stopAndGetElapsed(), m_isTask);
    if (m_isTask == false) {
      m_manager.popEvent(m_event);
    } else {
      m_manager.popContext(m_event);
    }
  }

  ScopeEvent(const ScopeEvent &) = delete;
  ScopeEvent &operator=(const ScopeEvent &) = delete;
  ScopeEvent(ScopeEvent &&) = delete;
  ScopeEvent &operator=(ScopeEvent &&) = delete;
};

#define ScopeEventUniqueKey_Core_To_Str_Ext(X) #X
#define ScopeEventUniqueKey_Core_To_Str(X)                                     \
  ScopeEventUniqueKey_Core_To_Str_Ext(X)
#define ScopeEventUniqueKey __FILE__ ScopeEventUniqueKey_Core_To_Str(__LINE__)

#define ScopeEventMultiRefKey std::string("-- multiref event --")
} // namespace util

#endif
