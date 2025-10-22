/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#pragma once

/** @file
 *  Implements binary output using MPI-IO.
 */

#include "ParticleRange.hpp"
#include "bonded_interactions/bonded_interaction_data.hpp"
#include "cell_system/CellStructure.hpp"

#include <string>

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

struct write_buffers {
  std::vector<double> pos;
  std::vector<double> vel;
  std::vector<int> id;
  std::vector<int> type;
};

/**
 * @brief Parallel binary output using MPI-IO.
 * To be called by all MPI processes. Aborts ESPResSo if an error occurs.
 * On 1 MPI rank, the error is converted to a runtime error and can be
 * recovered by removing any file that may have already been written.
 *
 * @param prefix Filepath prefix.
 * @param fields Specifier for which fields to dump.
 * @param bonded_ias Bonds to serialize.
 * @param particles Range of particles to serialize.
 * @param buffers Write buffers.
 */
void mpi_mpiio_common_write(std::string const &prefix, unsigned fields,
                            BondedInteractionsMap const &bonded_ias,
                            ParticleRange const &particles,
                            write_buffers &buffers);

/**
 * @brief Parallel binary input using MPI-IO.
 * To be called by all MPI processes. Aborts ESPResSo if an error occurs.
 * On 1 MPI rank, the error is converted to a runtime error and can be
 * recovered.
 *
 * @param prefix Filepath prefix.
 * @param fields Specifier for which fields to read.
 * @param cell_structure Bonds to serialize.
 */
void mpi_mpiio_common_read(std::string const &prefix, unsigned fields,
                           CellStructure &cell_structure);

} // namespace Mpiio
