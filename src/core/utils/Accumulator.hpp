#ifndef CORE_UTILS_ACCUMULATOR
#define CORE_UTILS_ACCUMULATOR

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/stream.hpp>

namespace Utils {

template <typename T> struct AccumulatorData {
  AccumulatorData()
      : m_mean_old(0.0), m_mean_new(0.0), m_variance_old(0.0),
        m_variance_new(0.0) {}
  T m_mean_old;
  T m_mean_new;
  T m_variance_old;
  T m_variance_new;
private:
  // Allow serialization to access non-public data members.
  friend class boost::serialization::access;

  template<typename Archive>
  void serialize(Archive& ar, const unsigned version) {
    ar & m_mean_old & m_mean_new & m_variance_old & m_variance_new;
  }
};

class Accumulator {
public:
  Accumulator(std::size_t N) : m_n(0), m_acc_data(N) {}
  void operator()(const std::vector<double>&);
  std::vector<double> get_mean() const;
  std::vector<double> get_variance() const;

private:
  std::size_t m_n;
  std::vector<AccumulatorData<double>> m_acc_data;
  // Allow serialization to access non-public data members.
  friend class boost::serialization::access;

  template<typename Archive>
  void serialize(Archive& ar, const unsigned version) {
    ar & m_n & m_acc_data;
  }
};

inline void Accumulator::operator()(const std::vector<double> &data) {
  if (data.size() != m_acc_data.size())
    throw std::runtime_error(
        "The given data size does not fit the initialized size!");
  ++m_n;
  if (m_n == 1) {
    auto acc_data_it = m_acc_data.begin();
    auto data_it = data.begin();
    for (; acc_data_it != m_acc_data.end() and data_it != data.end();
         ++acc_data_it, ++data_it) {
      (*acc_data_it).m_mean_old = (*acc_data_it).m_mean_new = *data_it;
      (*acc_data_it).m_variance_old = (*acc_data_it).m_variance_new = 0.0;
    }
  } else {
    // Calculate mean value.
    auto acc_data_it = m_acc_data.begin();
    auto data_it = data.begin();
    for (; acc_data_it != m_acc_data.end() and data_it != data.end();
         ++acc_data_it, ++data_it) {
      (*acc_data_it).m_mean_new = (*acc_data_it).m_mean_old +
                                  (*data_it - (*acc_data_it).m_mean_old) / m_n;
    }
    // Calculate variance.
    acc_data_it = m_acc_data.begin();
    data_it = data.begin();
    for (; acc_data_it != m_acc_data.end() and data_it != data.end();
         ++acc_data_it, ++data_it) {
      (*acc_data_it).m_variance_new =
          ((m_n - 1) * (*acc_data_it).m_variance_old +
           (*data_it - (*acc_data_it).m_mean_old) *
               (*data_it - (*acc_data_it).m_mean_new)) /
          m_n;
    }
    std::for_each(m_acc_data.begin(), m_acc_data.end(),
                   [](AccumulatorData<double> &acc_data) {
                     acc_data.m_mean_old = acc_data.m_mean_new;
                     acc_data.m_variance_old = acc_data.m_variance_new;
                   });
  }
}

inline std::vector<double> Accumulator::get_mean() const {
  std::vector<double> res;
  std::transform(
      m_acc_data.begin(), m_acc_data.end(), std::back_inserter(res),
      [](const AccumulatorData<double> &acc_data) { return acc_data.m_mean_new; });
  return res;
}

inline std::vector<double> Accumulator::get_variance() const {
  std::vector<double> res;
  std::transform(m_acc_data.begin(), m_acc_data.end(), std::back_inserter(res),
                 [](const AccumulatorData<double> &acc_data) {
                   return acc_data.m_variance_new;
                 });
  return res;
}


}

#endif
