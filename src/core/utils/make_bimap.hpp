/*
  Copyright (C) 2016 The ESPResSo project
  
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

#ifndef UTILS_MAKE_BIMAP_HPP
#define UTILS_MAKE_BIMAP_HPP

#include <boost/bimap.hpp>

namespace Utils {

/**
 * @brief Make a boost::bimap from a braced initializer list.
 *
 * Source: https://stackoverflow.com/a/31841462/3198615
 */
template <typename L, typename R>
boost::bimap<L, R>
make_bimap(std::initializer_list<typename boost::bimap<L, R>::value_type> list)
{
  return boost::bimap<L, R>(list.begin(), list.end());
}

} /* namespace Utils */

#endif
