/*
 * Copyright (C) 2016-2019 The ESPResSo project
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

#ifndef OBSERVABLES_PIDOBSERVABLE_HPP
#define OBSERVABLES_PIDOBSERVABLE_HPP

#include "Observable.hpp"
#include "Particle.hpp"

#include <utils/Span.hpp>

#include <vector>

namespace Observables {

/** %Particle-based observable.
 *
 *  Base class for observables extracting raw data from particle subsets and
 *  returning either the data or a statistic derived from it.
 */
class PidObservable : virtual public Observable {
  /** Identifiers of particles measured by this observable */
  std::vector<int> m_ids;

  virtual std::vector<double>
  evaluate(Utils::Span<const Particle *const> particles) const = 0;

public:
  explicit PidObservable(std::vector<int> ids) : m_ids(std::move(ids)) {}
  std::vector<double> operator()() const final;

  std::vector<int> &ids() { return m_ids; }
  std::vector<int> const &ids() const { return m_ids; }
};

} // Namespace Observables
#endif
