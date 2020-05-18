/*
 *
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
#ifndef NSQUARE_H
#define NSQUARE_H
/** \file
 *  This file contains the code for a simple n^2 particle loop.
 *
 *  The nsquare cell system performs a full n^2 particle interaction
 *  calculation over the simulation box. Therefore every node just
 *  has a single cell containing all local particles plus one ghost
 *  cell per other node. The communication is done via broadcasts
 *  (exchange_ghosts and update_ghosts) and reduce operations
 *  (collect_ghost_force).
 *
 *  The algorithm used for interaction calculation is parallelized,
 *  but a full communication is needed every time step. Let us assume
 *  that the number of nodes P is odd. Then a node p will do the
 *  interaction calculation with another node q iff \f$(q-p)\,mod\,P\f$ is even.
 *  Of course then every node has to do the same amount of work
 *  (provided the particles are distributed equally), and each
 *  interaction pair is done exactly once, since for odd P \f$r\,mod\,P\f$
 *  odd iff \f$-r\,mod\,P\f$ is even. For an even number of nodes,
 *  a virtual additional processor is assumed, with which no interactions occur.
 *  This means, that each communication cycle, 1 processor is idle, which is
 *  pretty ineffective.
 *
 *  The calculations themselves are just simple loops over all
 *  appropriate particle pairs.
 *
 *  Implementation in nsquare.cpp.
 */
#include "Cell.hpp"
#include "ParticleList.hpp"

#include <boost/mpi/communicator.hpp>

#include <vector>

/** setup the nsquare topology */
void nsq_topology_init(const boost::mpi::communicator &comm);

/** implements the load balancing as described above. */
void nsq_exchange_particles(int global_flag, ParticleList *displaced_parts,
                            std::vector<Cell *> &modified_cells);
#endif
