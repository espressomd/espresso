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
#ifndef OBSERVABLES_CYLINDRICALLBFLUXDENSITYPROFILEATPARTICLEPOSITIONS_HPP
#define OBSERVABLES_CYLINDRICALLBFLUXDENSITYPROFILEATPARTICLEPOSITIONS_HPP

#include "CylindricalPidProfileObservable.hpp"

#include "Particle.hpp"
#include "observables/ParticleTraits.hpp"

#include <utils/Span.hpp>

#include <cstddef>
#include <functional>
#include <vector>

namespace Observables {
class CylindricalLBFluxDensityProfileAtParticlePositions
    : public CylindricalPidProfileObservable {
public:
  using CylindricalPidProfileObservable::CylindricalPidProfileObservable;

  std::vector<double>
  evaluate(Utils::Span<std::reference_wrapper<const Particle>> particles,
           const ParticleObservables::traits<Particle> &traits) const override;

  std::vector<size_t> shape() const override {
    auto const b = n_bins();
    return {b[0], b[1], b[2], 3};
  }
};

} // Namespace Observables

#endif
