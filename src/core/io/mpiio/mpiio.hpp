/*
 * Copyright (C) 2010-2019 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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
#ifndef CORE_IO_MPIIO_MPIIO_HPP
#define CORE_IO_MPIIO_MPIIO_HPP

/** @file
 *  Implements binary output using MPI-IO.
 */

#include "ParticleRange.hpp"

namespace Mpiio {

/**
 * @brief Constants which indicate what to output.
 * To indicate the output of multiple fields, OR the
 * corresponding values.
 */
enum MPIIOOutputFields : unsigned int {
  MPIIO_OUT_NON = 0u,
  MPIIO_OUT_POS = 1u,
  MPIIO_OUT_VEL = 2u,
  MPIIO_OUT_TYP = 4u,
  MPIIO_OUT_BND = 8u,
};

/**
 * @brief Parallel binary output using MPI-IO.
 * To be called by all MPI processes. Aborts ESPResSo if an error occurs.
 * On 1 MPI rank, the error is converted to a runtime error and can be
 * recovered by removing any file that may have already been written.
 *
 * @param prefix Filepath prefix.
 * @param fields Specifier for which fields to dump.
 * @param particles Range of particles to serialize.
 */
void mpi_mpiio_common_write(const std::string &prefix, unsigned fields,
                            const ParticleRange &particles);

/**
 * @brief Parallel binary input using MPI-IO.
 * To be called by all MPI processes. Aborts ESPResSo if an error occurs.
 * On 1 MPI rank, the error is converted to a runtime error and can be
 * recovered.
 *
 * @param prefix Filepath prefix.
 * @param fields Specifier for which fields to read.
 */
void mpi_mpiio_common_read(const std::string &prefix, unsigned fields);

} // namespace Mpiio

#endif
