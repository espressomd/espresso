/*
 * Copyright (C) 2022-2023 The ESPResSo project
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <algorithm>
#include <cassert>
#include <memory>
#include <vector>

namespace EK {

template <class EKReaction> class EKReactions {
  using container_type = std::vector<std::shared_ptr<EKReaction>>;

public:
  using value_type = typename container_type::value_type;
  using iterator = typename container_type::iterator;
  using const_iterator = typename container_type::const_iterator;

private:
  container_type m_ekreactions;

public:
  bool contains(std::shared_ptr<EKReaction> const &ek_reaction) const noexcept {
    return std::find(begin(), end(), ek_reaction) != end();
  }
  void add(std::shared_ptr<EKReaction> const &ek_reaction) {
    assert(not contains(ek_reaction));
    m_ekreactions.emplace_back(ek_reaction);
  }
  void remove(std::shared_ptr<EKReaction> const &ek_reaction) {
    assert(contains(ek_reaction));
    m_ekreactions.erase(std::remove(begin(), end(), ek_reaction), end());
  }

  iterator begin() { return m_ekreactions.begin(); }
  iterator end() { return m_ekreactions.end(); }
  const_iterator begin() const { return m_ekreactions.begin(); }
  const_iterator end() const { return m_ekreactions.end(); }
  [[nodiscard]] bool empty() const { return m_ekreactions.empty(); }
};

} // namespace EK
