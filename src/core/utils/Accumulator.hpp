#ifndef CORE_UTILS_ACCUMULATOR
#define CORE_UTILS_ACCUMULATOR

namespace Utils {

class Accumulator {
public:
  Accumulator(std::size_t N)
      : m_n(0),
        m_mean_old(std::vector<double>(N, 0.0)),
        m_mean_new(std::vector<double>(N, 0.0)),
        m_variance_old(std::vector<double>(N, 0.0))
        m_variance_new(std::vector<double>(N, 0.0))
        {}
  void operator()(std::vector<double>);
  std::vector<double> get_mean() const;
  std::vector<double> get_variance() const;

private:
  std::size_t m_n;
  std::vector<double> m_mean_old;
  std::vector<double> m_mean_new;
  std::vector<double> m_variance_old;
  std::vector<double> m_variance_new;
}
}

#endif
