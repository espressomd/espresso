/*
Copyright (C) 2010-2018 The ESPResSo project

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
#ifndef CORE_UTILS_SERIALIZATION_IBM_CUDA_PARTICLE_DATA_INPUT_HPP
#define CORE_UTILS_SERIALIZATION_IBM_CUDA_PARTICLE_DATA_INPUT_HPP

#include "virtual_sites/lb_inertialess_tracers_cuda_interface.hpp"

namespace boost {
namespace serialization {
template <class Archive>
void serialize(Archive &ar, IBM_CUDA_ParticleDataInput &d,
               const unsigned int /* version */) {
  ar &d.pos;
  ar &d.f;
  ar &d.is_virtual;
}

template <class Archive>
void serialize(Archive &ar, IBM_CUDA_ParticleDataOutput &d,
               const unsigned int /* version */) {
  ar &d.v;
}
} // namespace serialization
} // namespace boost

#endif
