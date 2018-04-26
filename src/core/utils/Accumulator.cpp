#include "Accumulator.hpp"

namespace Utils {

Accumulator::operator()(std::vector<double> data) {
  if (data.size() != m_mean_old.size())
    throw std::runtime_error(
        "The given data size does not fit the initialized size!");
  ++m_n;
  if (m_n == 1) {
    m_mean_old = m_mean_new = data;
  } else {
    std::transform(m_mean_old.begin(), m_mean_old.end(), data.begin(),
                   m_mean_new.begin(),
                   [m_n](double const &mean_old, double const &datasample) {
                     return (mean_old + (datasample - mean_old) / m_n)
                   });
    for (auto variance_new_it = m_variance_new.begin(),
              auto variance_old_it = m_variance_old.cbegin(),
              auto datasample_it = data.cbegin(),
              auto mean_new_it = m_mean_new.cbegin(),
              auto mean_old_it = m_mean_old.cbegin();
         variance_new_it != m_variance_new.end() and
         variance_old_it != m_variance_old.end() and
         datasample_it != data.end() and mean_new_it != m_mean_new.end() and
         mean_old_it != m_mean_old.end();
         ++variance_new_it, ++variance_old_it, ++datasample_it, ++mean_new_it,
              ++mean_old_it) {
      *variance_new_it =
          ((m_n - 1) * *variance_old_it +
           (*datasample_it - *mean_old_it) * (*datasample_it - *mean_new_it)) /
          m_n;
    }
    m_mean_old = m_mean_new;
    m_variance_old = m_mean_new;
  }
}
