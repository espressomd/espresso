/*
 * Copyright (C) 2016-2022 The ESPResSo project
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

#include "Observable.hpp"

#include "PidObservable.hpp"

#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace Observables {

/** Radial distribution function.
 */
class RDF : public Observable {
  /** Identifiers of the reference particles */
  std::vector<int> m_ids1;
  /** Identifiers of the distant particles */
  std::vector<int> m_ids2;

  std::vector<double>
  evaluate(boost::mpi::communicator const &comm,
           ParticleReferenceRange const &local_particles_1,
           ParticleReferenceRange const &local_particles_2,
           const ParticleObservables::traits<Particle> &traits) const;

public:
  // Range of the profile.
  double min_r, max_r;
  // Number of bins
  std::size_t n_r_bins;

  std::vector<std::size_t> shape() const override { return {n_r_bins}; }

  RDF(std::vector<int> ids1, std::vector<int> ids2, int n_r_bins, double min_r,
      double max_r)
      : m_ids1(std::move(ids1)), m_ids2(std::move(ids2)), min_r(min_r),
        max_r(max_r) {
    if (max_r <= min_r)
      throw std::runtime_error("max_r has to be > min_r");
    if (n_r_bins <= 0)
      throw std::domain_error("n_r_bins has to be >= 1");
    this->n_r_bins = static_cast<std::size_t>(n_r_bins);
  }
  std::vector<double>
  operator()(boost::mpi::communicator const &comm) const final;

  std::vector<int> &ids1() { return m_ids1; }
  std::vector<int> &ids2() { return m_ids2; }
  std::vector<int> const &ids1() const { return m_ids1; }
  std::vector<int> const &ids2() const { return m_ids2; }
};

} // Namespace Observables
