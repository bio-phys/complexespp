///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef INATIMER_HPP
#define INATIMER_HPP

#include <chrono>


class InaTimer {
    using double_second_time = std::chrono::duration< double, std::ratio< 1, 1 > >;

    std::chrono::high_resolution_clock::time_point m_start; ///< m_start time (start)
    std::chrono::high_resolution_clock::time_point m_end;   ///< stop time (stop)
    std::chrono::nanoseconds m_cumulate;                    ///< the m_cumulate time

public:
    /// Constructor
    InaTimer() {
        start();
    }

    /// Copy constructor
    InaTimer(const InaTimer& other) = delete;
    /// Copies an other timer
    InaTimer& operator=(const InaTimer& other) = delete;
    /// Move constructor
    InaTimer(InaTimer&& other) = delete;
    /// Copies an other timer
    InaTimer& operator=(InaTimer&& other) = delete;

    /** Rest all the values, and apply start */
    void reset() {
        m_start    = std::chrono::high_resolution_clock::time_point();
        m_end      = std::chrono::high_resolution_clock::time_point();
        m_cumulate = std::chrono::nanoseconds();
        start();
    }

    /** Start the timer */
    void start() {
        m_start = std::chrono::high_resolution_clock::now();
    }

    /** Stop the current timer */
    void stop() {
        m_end = std::chrono::high_resolution_clock::now();
        m_cumulate +=
            std::chrono::duration_cast< std::chrono::nanoseconds >(m_end - m_start);
    }

    /** Return the elapsed time between start and stop (in second) */
    double getElapsed() const {
        return std::chrono::duration_cast< double_second_time >(
                   std::chrono::duration_cast< std::chrono::nanoseconds >(
                       m_end - m_start))
            .count();
    }

    /** Return the total counted time */
    double getCumulated() const {
        return std::chrono::duration_cast< double_second_time >(m_cumulate).count();
    }

    /** End the current counter (stop) and return the elapsed time */
    double stopAndGetElapsed() {
        stop();
        return getElapsed();
    }
};

#endif
