/*
 * Copyright (C) 2010-2020 The ESPResSo project
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
#ifndef SERIALIZATION_SD_PARTICLE_DATA_HPP
#define SERIALIZATION_SD_PARTICLE_DATA_HPP

#include "stokesian_dynamics/sd_interface.hpp"

#include <boost/serialization/is_bitwise_serializable.hpp>
#include <boost/serialization/split_free.hpp>

#if defined(STOKESIAN_DYNAMICS) or defined(STOKESIAN_DYNAMICS_GPU)

BOOST_IS_BITWISE_SERIALIZABLE(SD_particle_data)

namespace boost {
namespace serialization {
template <typename Archive>
void load(Archive &ar, SD_particle_data &p,
          const unsigned int /* file_version */) {
  ar >> make_array(reinterpret_cast<char *>(&p), sizeof(SD_particle_data));
}

template <typename Archive>
void save(Archive &ar, SD_particle_data const &p,
          const unsigned int /* file_version */) {
  /* Cruel but effective */
  ar << make_array(reinterpret_cast<char const *>(&p),
                   sizeof(SD_particle_data));
}

template <class Archive>
void serialize(Archive &ar, SD_particle_data &p,
               const unsigned int file_version) {
  split_free(ar, p, file_version);
}
} // namespace serialization
} // namespace boost

#endif
#endif
