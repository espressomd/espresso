/*
  Copyright (C) 2010,2012,2013,2014,2015,2016 The ESPResSo project
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

#include "RunningAverage.hpp"

namespace Utils {
namespace Statistics {

template<typename Scalar>
void RunningAverage<Scalar>::add_sample(Scalar s) {
  m_n++;

  if(m_n == 1) {
    m_old_avg = m_new_avg = s;
    m_old_var = 0.0;
  } else {
    m_new_avg = m_old_avg + (s - m_old_avg)/m_n;
    m_new_var = m_old_var + (s - m_old_avg)*(s - m_new_avg);

    m_old_avg = m_new_avg;
    m_old_var = m_new_var;
  }
}

template<typename Scalar>
Scalar RunningAverage<Scalar>::avg() const { 
    if(m_n > 0)
      return m_new_avg;
    else
      return 0.0;
  }

template<typename Scalar>
Scalar RunningAverage<Scalar>::var() const {
  if(m_n > 1)
    return m_new_var / (m_n - 1.);
  else
    return 0.0;
}

}
}
