/*
  Copyright (C) 2010-2018 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
    Max-Planck-Institute for Polymer Research, Theory Group

  This file is part of ESPResSo.

  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef __RUNING_AVERAGE_HPP
#define __RUNING_AVERAGE_HPP

#include <algorithm>
#include <cmath>
#include <limits>

namespace Utils {
namespace Statistics {

/**
 * \brief Keep running average and variance.
 * The average should be numerically stable.
 */
template <typename Scalar> class RunningAverage {
public:
  RunningAverage()
      : m_n(0), m_old_avg(0), m_new_avg(0), m_old_var(0), m_new_var(0.0),
        m_min(std::numeric_limits<Scalar>::max()),
        m_max(-std::numeric_limits<Scalar>::max()) {}
  void add_sample(Scalar s) {
    m_n++;

    if (m_n == 1) {
      m_old_avg = m_new_avg = s;
    } else {
      m_new_avg = m_old_avg + (s - m_old_avg) / m_n;
      m_new_var = m_old_var + (s - m_old_avg) * (s - m_new_avg);

      m_old_avg = m_new_avg;
      m_old_var = m_new_var;
    }

    m_min = std::min(m_min, s);
    m_max = std::max(m_max, s);
  }

  void clear() { m_n = 0; }

  int n() const { return m_n; }

  /** Average of the samples */
  Scalar avg() const {
    if (m_n > 0)
      return m_new_avg;

    return 0.0;
  }
  /** Variance of the samples */
  Scalar var() const {
    if (m_n > 1)
      return m_new_var / m_n;

    return 0.0;
  }

  /** Standard deviation of the samples */
  Scalar sig() const { return std::sqrt(var()); }

  /** Minimum */
  Scalar min() const { return m_min; }

  /** Minimum */
  Scalar max() const { return m_max; }

private:
  int m_n;
  Scalar m_old_avg, m_new_avg;
  Scalar m_old_var, m_new_var;
  Scalar m_min, m_max;
};
} // namespace Statistics
} // namespace Utils

#endif
