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

#ifndef _GENERIC_DD_HPP
#define _GENERIC_DD_HPP

/** @file generic_dd.hpp
 * This file implements a more generic way of ESPResSo's own "domain
 * decomposition" cell system. "Generic" means it allows for non-Cartesian
 * decompositions of the domain. The underlying parallel linked cell grid
 * is deferred to "librepa", a library for parallel grid implementations. These
 * functions only adapt ESPResSo's internal datastructures and provide
 * functionality that ESPResSo expects based on the information from librepa.
 */
#include <string>
#include <vector>

#include "cells.hpp"
#include "metric.hpp"

namespace generic_dd {

/** @brief Adjust the domain decomposition to a change in the geometry.
 *
 * May call cells_re_init if the new range cannot be accounted for with the
 * existing grid.
 *
 * @see dd_on_geometry_change
 * @param flags flags for the geometry change
 * @param range Minimum required cell width
 */
void on_geometry_change(int flags, double range);

/** Initialize the generic-dd cell system.
 * After calling this function, all other functions exported from the generic-dd
 * module can safely be called.
 *
 * @see topology_release
 * @param range Minimum required cell width.
 * @param is_repart If false triggers a creation of a new grid.
 */
void topology_init(double range, bool is_repart = false);

/** Frees data allocated by the generic-dd cell system.
 * After calling this, do not call any other function exported from the
 * generic-dd module *but* topology_init.
 */
void topology_release();

/** Exchanges particles and tracks the cells it changed.
 *
 * Particles not belonging to this node are communicated to the corresponding
 * new owners. If global_flag is false, this happens only with
 * the processes neighboring the current one.
 *
 * @param global_flag Indicates global or local communication
 * @param displaced_particles A list of particles to be communicated.
 * @param modified_cells Local cells that are modified in the course of this
 * function are pushed back here.
 */
void exchange_particles(int global_flag, ParticleList *displaced_particles,
                        std::vector<Cell *> &modified_cells);

/** Use grid specified by "desc" as underlying parallel grid implementation in
 * subsequent calls to topology_init.
 */
void set_grid(const std::string &desc);

/** Repartitions the grid unconditionally if the underlying grid implementation
 * does support repartitioning.
 * After this function returns, the internal state of ESPResSo is valid, again,
 * and the simulation may continue (with the new subdomain layout).
 *
 * @param m Metric to use.
 */
void repartition(const repart::Metric &m);

/** Deliver an implementation-defined command to the current partitioner/grid.
 */
void command(const std::string &cmd);

/** Returns a vector of supported grid types by librepa.
 */
std::vector<std::string> librepa_supported_grid_types();

} // namespace generic_dd

#endif
