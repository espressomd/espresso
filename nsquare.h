// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2004; all rights reserved unless otherwise stated.
#ifndef NSQUARE_H
#define NSQUARE_H
/** \file nsquare.h
    This file contains the code for a simple n^2 particle loop.

    <b>Responsible:</b>
    <a href="mailto:arnolda@mpip-mainz.mpg.de">Axel</a>

    The nsquare cell system performs a full n^2 particle interaction
    calculation over the simulation box.  Therefore every node just
    has a single cell containing all local particles plus one ghost
    cell per other node. The communication is done via broadcasts
    (exchange_ghosts and update_ghosts) and reduce operations
    (collect_ghost_force).
    
    The algorithm used for interaction calculation is parallelized,
    but a full communication is needed every time step. Let us assume
    that the number of nodes P is odd. Then a node p will do the
    interaction calculation with another node q iff \f$(q-p)\,mod\,P\f$ is even.
    Of course then every node has to do the same amount of work
    (provided the particles are distributed equally), and each
    interaction pair is done exactly once, since for odd P \f$r\,mod\,P\f$
    odd iff \f$-r\,mod\,P\f$ is even. For an even number of nodes,
    a virtual additional processor is assumed, with which no interactions occur.
    This means, that each communication cycle, 1 processor is idle, which is
    pretty ineffective.

    The second main part of this cell system is a load balancer which
    at the beginning of the integration is called and balances the
    numbers of particles between the nodes. This means that the number
    of particles per node should be around \f$N/P\f$, so the goal is to
    let every processor have a number of particles which is one of the
    two integers closest to \f$N/P\f$. The algorithm is greedy, i. e. it
    searches for the node with most and least particles and transfers
    as much as possible particles between them so that at least one of
    them satisfies the constraints. Of course the algorithm terminates
    if both satisfy the condition without transfer.

    The calculations themselves are just simple loops over all
    appropriate particle pairs.
*/

#include "cells.h"

/** always returns the one local cell */
Cell *nsq_position_to_cell(double pos[3]);

/** always returns the one local cell */
void nsq_topology_release();

/** setup the nsquare topology */
void nsq_topology_init(CellPList *local);

/** implements the load balancing as described above. */
void nsq_balance_particles();

/** n^2 force calculation */
void nsq_calculate_ia();

/** n^2 energy calculation */
void nsq_calculate_energies();

/** n^2 pressure calculation */
void nsq_calculate_virials();

#endif
