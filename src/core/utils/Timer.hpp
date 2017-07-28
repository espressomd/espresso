#ifndef UTILS_TIMER_HPP
#define UTILS_TIMER_HPP

#include <map>
#include <string>
#include <unordered_map>

#include <boost/serialization/access.hpp>

#include <mpi.h>

#include "MpiCallbacks.hpp"
#include "utils/statistics/RunningAverage.hpp"

/**
 * @brief Time events and keep staistics.
 *
 * Times are wall time and in seconds.
 */
namespace Utils {
namespace Timing {
class Timer {
public:
  struct Stats {
    Stats() {}
    Stats(double avg, double sig, double var, double min, double max, int n)
        : m_avg(avg), m_sig(sig), m_var(var), m_min(min), m_max(max), m_n(n) {}

    double avg() const { return m_avg; }

    double sig() const { return m_sig; }

    double var() const { return m_var; }

    double min() const { return m_min; }

    double max() const { return m_max; }

    int n() const { return m_n; }

  private:
    double m_avg, m_sig, m_var, m_min, m_max;
    int m_n;

    friend boost::serialization::access;
    template <class Archive>
    void serialize(Archive &ar, const unsigned int /* version */) {
      ar &m_avg;
      ar &m_sig;
      ar &m_var;
      ar &m_min;
      ar &m_max;
      ar &m_n;
    }
  };

  void start() { m_mark = MPI_Wtime(); }

  double stop() {
    const double time = MPI_Wtime() - m_mark;
    if (time > 0.0)
      m_running_average.add_sample(time);

    return time;
  }

  Stats stats() const {
    return Stats(m_running_average.avg(), m_running_average.sig(),
                 m_running_average.var(), m_running_average.min(),
                 m_running_average.max(), m_running_average.n());
  }

  void reset() { m_running_average.clear(); }

  static Timer &get_timer(const std::string &name) { return m_timers[name]; }

  static std::map<std::string, Stats> get_stats() {
    std::map<std::string, Stats> ret;

    for (auto &it : m_timers) {
      ret[it.first] = it.second.stats();
    }

    return ret;
  }

  /**
   * @brief Collect stats from all nodes.
   *
   * Can only be called on the master node.
   */
  static std::vector<std::map<std::string, Utils::Timing::Timer::Stats>>
    get_all_stats();
private:
  friend Communication::MpiCallbacks;
  static void mpi_callback(int, int);

  static const Communication::CallbackAdder cb_adder;

  static std::unordered_map<std::string, Timer> m_timers;
  Statistics::RunningAverage<double> m_running_average;
  double m_mark;
};
}
}

#endif
