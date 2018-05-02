#ifndef CORE_UTILS_ACCUMULATOR
#define CORE_UTILS_ACCUMULATOR

#include <boost/serialization/access.hpp>

namespace Utils {

template <typename T> struct AccumulatorData {
  AccumulatorData() = default;
  T mean;
  T variance;

private:
  // Allow serialization to access non-public data members.
  friend class boost::serialization::access;

  template <typename Archive>
  void serialize(Archive &ar, const unsigned version) {
    ar &mean &variance;
  }
};

class Accumulator {
public:
  Accumulator(std::size_t N) : m_n(0), m_acc_data(N) {}
  void operator()(const std::vector<double> &);
  std::vector<double> get_mean() const;
  std::vector<double> get_variance() const;

private:
  std::size_t m_n;
  std::vector<AccumulatorData<double>> m_acc_data;
  // Allow serialization to access non-public data members.
  friend class boost::serialization::access;

  template <typename Archive>
  void serialize(Archive &ar, const unsigned version) {
    ar &m_n &m_acc_data;
  }
};

inline void Accumulator::operator()(const std::vector<double> &data) {
  if (data.size() != m_acc_data.size())
    throw std::runtime_error(
        "The given data size does not fit the initialized size!");
  ++m_n;
  if (m_n == 1) {
    std::transform(data.begin(), data.end(), m_acc_data.begin(),
                   [](double d) -> AccumulatorData<double> {
                     return {d, 0.0};
                   });
  } else {
    std::transform(
        m_acc_data.begin(), m_acc_data.end(), data.begin(), m_acc_data.begin(),
        [this](AccumulatorData<double> &a,
               double d) -> AccumulatorData<double> {
          auto const old_mean = a.mean;
          auto const new_mean = old_mean + (d - old_mean) / m_n;
          auto const new_variance =
              ((m_n - 1) * a.variance + (d - old_mean) * (d - new_mean)) / m_n;
          return {new_mean, new_variance};
        });
  }
}

inline std::vector<double> Accumulator::get_mean() const {
  std::vector<double> res;
  std::transform(
      m_acc_data.begin(), m_acc_data.end(), std::back_inserter(res),
      [](const AccumulatorData<double> &acc_data) { return acc_data.mean; });
  return res;
}

inline std::vector<double> Accumulator::get_variance() const {
  std::vector<double> res;
  std::transform(m_acc_data.begin(), m_acc_data.end(), std::back_inserter(res),
                 [](const AccumulatorData<double> &acc_data) {
                   return acc_data.variance;
                 });
  return res;
}
}

#endif
