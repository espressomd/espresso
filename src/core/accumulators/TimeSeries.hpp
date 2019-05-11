//
// Created by florian on 11.05.19.
//

#ifndef ESPRESSO_TIMESERIES_HPP
#define ESPRESSO_TIMESERIES_HPP

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
  // The accumulator struct has to be initialized with the correct vector size,
  // therefore the order of init is important.
  TimeSeries(std::shared_ptr<Observables::Observable> const &obs, int delta_N)
      : AccumulatorBase(delta_N), m_obs(obs) {}

  void update() override;
  /* Partial serialization of state that is not accessible
     via the interface. */
  std::string get_internal_state() const;
  void set_internal_state(std::string const &);

  const std::vector<std::vector<double>> &time_series() const { return m_data; }
  void clear() {
    m_data.clear();
  }

private:
  std::shared_ptr<Observables::Observable> m_obs;
  std::vector<std::vector<double>> m_data;
};

} // namespace Accumulators

#endif // ESPRESSO_TIMESERIES_HPP
