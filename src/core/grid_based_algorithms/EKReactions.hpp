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

#ifndef ESPRESSO_EKREACTIONS_HPP
#define ESPRESSO_EKREACTIONS_HPP

#include <algorithm>
#include <memory>
#include <vector>

template <class EKReaction> class EKReactions {
  using container_type = std::vector<std::shared_ptr<EKReaction>>;

public:
  using value_type = typename container_type::value_type;
  using iterator = typename container_type::iterator;
  using const_iterator = typename container_type::const_iterator;

private:
  container_type m_ekreactions;

public:
  void add(std::shared_ptr<EKReaction> const &c) {
    assert(std::find(m_ekreactions.begin(), m_ekreactions.end(), c) ==
           m_ekreactions.end());

    m_ekreactions.emplace_back(c);
  }
  void remove(std::shared_ptr<EKReaction> const &c) {
    assert(std::find(m_ekreactions.begin(), m_ekreactions.end(), c) !=
           m_ekreactions.end());
    m_ekreactions.erase(
        std::remove(m_ekreactions.begin(), m_ekreactions.end(), c),
        m_ekreactions.end());
  }

  iterator begin() { return m_ekreactions.begin(); }
  iterator end() { return m_ekreactions.end(); }
  const_iterator begin() const { return m_ekreactions.begin(); }
  const_iterator end() const { return m_ekreactions.end(); }
  [[nodiscard]] bool empty() const { return m_ekreactions.empty(); }
};

#endif
