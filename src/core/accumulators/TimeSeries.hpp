#ifndef CORE_ACCUMULATORS_TIMESERIES_HPP
#define CORE_ACCUMULATORS_TIMESERIES_HPP

#include "AccumulatorBase.hpp"
#include "observables/Observable.hpp"

namespace Accumulators {

/**
 * @brief Record values of an observable.
 *
 * This is a very simple accumulator that stores
 * the current value of an observable everytime
 * it is updated.
 *
 */
class TimeSeries : public AccumulatorBase {
public:
  TimeSeries(std::shared_ptr<Observables::Observable> obs, int delta_N)
      : AccumulatorBase(delta_N), m_obs(std::move(obs)) {}

  void update() override;
  std::string get_internal_state() const;
  void set_internal_state(std::string const &);

  const std::vector<std::vector<double>> &time_series() const { return m_data; }
  void clear() { m_data.clear(); }

private:
  std::shared_ptr<Observables::Observable> m_obs;
  std::vector<std::vector<double>> m_data;
};

} // namespace Accumulators

#endif
