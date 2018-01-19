/*
  Copyright (C) 2017 The ESPResSo project

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

#ifndef UTILS_KEYS_HPP
#define UTILS_KEYS_HPP

#include <vector>
#include <algorithm>

namespace Utils {
  /**
   * @brief Return the keys of a map type.
   *
   * Returns a vector of copies of the keys
   * of a map, unordered_map, ...
   */
  template <typename Map>
  auto keys(Map const &m) -> std::vector<typename Map::key_type> {
    using value_type = typename Map::value_type;
    using std::begin;
    using std::end;

    std::vector<typename Map::key_type> ret(m.size());

    std::transform(begin(m), end(m), ret.begin(),
                   [](value_type const &kv) { return kv.first; });
    return ret;
  }
}

#endif
