#include "TimeSeries.hpp"

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/serialization/vector.hpp>

namespace Accumulators {
void TimeSeries::update() { m_data.emplace_back(m_obs->operator()()); }

std::string TimeSeries::get_internal_state() const {
  std::stringstream ss;
  boost::archive::binary_oarchive oa(ss);

  oa << m_data;

  return ss.str();
}

void TimeSeries::set_internal_state(std::string const &state) {
  namespace iostreams = boost::iostreams;
  iostreams::array_source src(state.data(), state.size());
  iostreams::stream<iostreams::array_source> ss(src);
  boost::archive::binary_iarchive ia(ss);

  ia >> m_data;
}
} // namespace Accumulators
