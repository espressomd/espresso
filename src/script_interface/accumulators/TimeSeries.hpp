#ifndef SCRIPT_INTERFACE_ACCUMULATORS_TIME_SERIES_HPP
#define SCRIPT_INTERFACE_ACCUMULATORS_TIME_SERIES_HPP

#include "AccumulatorBase.hpp"
#include "core/accumulators/TimeSeries.hpp"
#include "observables/Observable.hpp"

#include <boost/range/algorithm/transform.hpp>
#include <utils/as_const.hpp>

#include <memory>

namespace ScriptInterface {
namespace Accumulators {

class TimeSeries : public AccumulatorBase {
public:
  /* as_const is to make obs read-only. */
  TimeSeries() { add_parameters({{"obs", Utils::as_const(m_obs)}}); }

  void construct(VariantMap const &params) override {
    set_from_args(m_obs, params, "obs");

    if (m_obs)
      m_accumulator = std::make_shared<::Accumulators::TimeSeries>(
          m_obs->observable(), get_value_or<int>(params, "delta_N", 1));
  }

  Variant call_method(std::string const &method,
                      VariantMap const &parameters) override {
    if (method == "update") {
      m_accumulator->update();
    }
    if (method == "time_series") {
      auto const &series = m_accumulator->time_series();
      std::vector<Variant> ret(series.size());

      boost::transform(
          series, ret.begin(),
          [](std::vector<double> const &sample) { return sample; });

      return ret;
    }
    if (method == "clear") {
      m_accumulator->clear();
    }

    return {};
  }

  Variant get_state() const override {
    std::vector<Variant> state(2);
    state[0] = ScriptInterfaceBase::get_state();
    state[1] = m_accumulator->get_internal_state();
    return state;
  }

  std::shared_ptr<::Accumulators::AccumulatorBase> accumulator() override {
    return std::static_pointer_cast<::Accumulators::AccumulatorBase>(
        m_accumulator);
  }

  std::shared_ptr<const ::Accumulators::AccumulatorBase>
  accumulator() const override {
    return std::static_pointer_cast<::Accumulators::AccumulatorBase>(
        m_accumulator);
  }

private:
  void set_state(Variant const &state) override {
    auto const &state_vec = boost::get<std::vector<Variant>>(state);
    ScriptInterfaceBase::set_state(state_vec.at(0));
    m_accumulator->set_internal_state(boost::get<std::string>(state_vec.at(1)));
  }

  /* The actual accumulator */
  std::shared_ptr<::Accumulators::TimeSeries> m_accumulator;
  std::shared_ptr<Observables::Observable> m_obs;
};

} // namespace Accumulators
} /* namespace ScriptInterface */

#endif // SCRIPT_INTERFACE_ACCUMULATORS_TIMESERIES_HPP
