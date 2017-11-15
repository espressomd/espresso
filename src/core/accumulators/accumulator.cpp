#include "accumulator.hpp"

namespace Accumulators {

void Accumulator::add_sample() {
  m_acc(std::vector(m_obs.last_value));
}

double Accumulator::get_mean() {
  boost::accumulators::mean(m_acc);
}

double Accumulator::get_variance() {
  boost::accumulators::variance(m_acc);
}
}
