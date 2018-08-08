/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
    Max-Planck-Institute for Polymer Research, Theory Group

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
/** \file reaction.cpp
 *
 */

#include "swimmer_reaction.hpp"
#include "cells.hpp"
#include "errorhandling.hpp"
#include "forces.hpp"
#include "initialize.hpp"
#include "utils.hpp"
#include "random.hpp"
#include "communication.hpp"
#include "integrate.hpp"
#include "grid.hpp"

#include <algorithm>
#include <random>
#include <vector>

reaction_struct reaction;

#ifdef SWIMMER_REACTIONS

void reactions_sanity_checks() {
  if (reaction.ct_rate != 0.0) {

    if (!cell_structure.use_verlet_list || cell_structure.type != CELL_STRUCTURE_DOMDEC) {
      runtimeErrorMsg() << "The SWIMMER_REACTIONS feature requires verlet "
                           "lists and domain decomposition";
    }

    if (max_cut < reaction.range) {
      runtimeErrorMsg() << "Reaction range of " << reaction.range
                        << " exceeds maximum cutoff of " << max_cut;
    }
  }
}

void local_setup_reaction() {

  /* Make available the various reaction parameters */
  MPI_Bcast(&reaction.reactant_type, 1, MPI_INT, 0, comm_cart);
  MPI_Bcast(&reaction.product_type, 1, MPI_INT, 0, comm_cart);
  MPI_Bcast(&reaction.catalyzer_type, 1, MPI_INT, 0, comm_cart);
  MPI_Bcast(&reaction.range, 1, MPI_DOUBLE, 0, comm_cart);
  MPI_Bcast(&reaction.ct_rate, 1, MPI_DOUBLE, 0, comm_cart);
  MPI_Bcast(&reaction.eq_rate, 1, MPI_DOUBLE, 0, comm_cart);
  MPI_Bcast(&reaction.sing_mult, 1, MPI_INT, 0, comm_cart);
  MPI_Bcast(&reaction.swap, 1, MPI_INT, 0, comm_cart);

  /* Create the various reaction related types (categories) */
  make_particle_type_exist(reaction.catalyzer_type);
  make_particle_type_exist(reaction.product_type);
  make_particle_type_exist(reaction.reactant_type);

  /* Make ESPResSo aware that reactants and catalyst are interacting species */
  IA_parameters *data =
      get_ia_param_safe(reaction.reactant_type, reaction.catalyzer_type);

  /* Used for the range of the verlet lists */
  data->REACTION_range = reaction.range;

  /* Broadcast interaction parameters */
  mpi_bcast_ia_params(reaction.reactant_type, reaction.catalyzer_type);
}

void integrate_reaction_noswap() {
  double dist2, vec21[3], ct_ratexp, eq_ratexp, rand, bernoulli;

  if (reaction.ct_rate > 0.0) {

    /* Determine the reaction rate */
    ct_ratexp = exp(-time_step * reaction.ct_rate);

    on_observable_calc();

    for (int c = 0; c < local_cells.n; c++) {

      /* Take into account only those cell neighbourhoods for which
         the central cell contains a catalyzer particle */

      auto check_catalyzer = 0;

      auto cell = local_cells.cell[c];
      auto p1 = cell->part;
      auto np = cell->n;

      for (int i = 0; i < np; i++) {
        if (p1[i].p.type == reaction.catalyzer_type) {
          check_catalyzer = 1;
          break;
        }
      }

      /* If the central cell contains a catalyzer particle, ...*/
      if (check_catalyzer != 0) {

        for (auto &pair : cell->m_verlet_list) {
          auto p1 = pair.first;      // pointer to particle 1
          auto p2 = pair.second; // pointer to particle 2

          if ((p1->p.type == reaction.reactant_type &&
               p2->p.type == reaction.catalyzer_type) ||
              (p2->p.type == reaction.reactant_type &&
               p1->p.type == reaction.catalyzer_type)) {
            get_mi_vector(vec21, p1->r.p, p2->r.p);
            dist2 = sqrlen(vec21);

            /* Count the number of times a reactant particle is
               checked against a catalyst */
            if (dist2 < reaction.range * reaction.range) {

              if (p1->p.type == reaction.reactant_type) {
                p1->p.catalyzer_count++;
              } else {
                p2->p.catalyzer_count++;
              }
            }
          }
        }
      }
    }

  /* Now carry out the reaction on the particles which are tagged */
  for (int c = 0; c < local_cells.n; c++) {
    auto cell = local_cells.cell[c];
    auto p1 = cell->part;
    auto np = cell->n;

    for (int i = 0; i < np; i++) {
      if (p1[i].p.type == reaction.reactant_type) {

        if (p1[i].p.catalyzer_count > 0) {

          if (reaction.sing_mult == 0) {
            rand = d_random();

            bernoulli = pow(ct_ratexp, p1[i].p.catalyzer_count);

            if (rand > bernoulli) {
              p1[i].p.type = reaction.product_type;
            }

          } else /* We only consider each reactant once */
          {
            rand = d_random();

            if (rand > ct_ratexp) {
              p1[i].p.type = reaction.product_type;
            }
          }

          p1[i].p.catalyzer_count = 0;
        }
      }
    }
    }

    /* We only need to do something when the equilibrium
       reaction rate constant is nonzero */
    if (reaction.eq_rate > 0.0) {

      eq_ratexp = exp(-time_step * reaction.eq_rate);

      for (int c = 0; c < local_cells.n; c++) {
        auto cell = local_cells.cell[c];
        auto p1 = cell->part;
        auto np = cell->n;

        /* Convert products into reactants and vice versa
           according to the specified rate constant */
        for (int i = 0; i < np; i++) {

          if (p1[i].p.type == reaction.product_type) {
            rand = d_random();

            if (rand > eq_ratexp) {
              p1[i].p.type = reaction.reactant_type;
            }
          } else if (p1[i].p.type == reaction.reactant_type) {
            rand = d_random();

            if (rand > eq_ratexp) {
              p1[i].p.type = reaction.product_type;
            }
          }
        }
      }
    }

    on_particle_change();
  }
}

#ifdef ROTATION

bool in_lower_half_space(Particle p1, Particle p2) {
  // This function determines whether the particle p2 is in the lower
  // half space of particle p1
  auto const distvec = get_mi_vector(p1.r.p, p2.r.p);
  double dot = p1.r.quatu * distvec;
  int sgn = Utils::sgn(dot);
  return (sgn + 1) / 2;
}

void integrate_reaction_swap() {
  int np;
  Particle *p_local, *p_neigh;
  Cell *cell;
  double dist2, vec21[3], ct_ratexp, rand;
  int n_reactions;

#ifdef ELECTROSTATICS
  double product_q = 0.0, reactant_q = 0.0;
#endif
  std::vector<int> catalyzers, reactants, products;

  // To randomize the lists of cells below we need a random number
  // engine.  In previous versions of C++ this was a single function
  // std::random_shuffle, which is unfortunately deprecated in C++17
  // and discouraged from C++11 on.
  std::random_device rd;  // non-deterministic random number generator
                          // using hardware entropy source
  std::mt19937 rng(rd()); // Mersenne Twister by Matsumoto and Nishimura

  // If multiple catalyzers get close to each other, they might eat up
  // each others reactants.  If we traverse the cells in a sorted
  // manner, then catalyzers in the upper left will use up all the
  // reactants of the catalyzers which are below right of them.  This
  // process is biased.  To rectify this issue we set up a vector
  // which goes through the cells in a randomized manner.
  std::vector<int> rand_cells(local_cells.n);
  for (int i = 0; i < local_cells.n; i++)
    rand_cells[i] = i;
  std::shuffle(rand_cells.begin(), rand_cells.end(), rng);

  if (reaction.ct_rate > 0.0) {
    // Determine the reaction rate
    ct_ratexp = exp(-time_step * reaction.ct_rate);

    on_observable_calc();

    // Iterate over all the local cells
    for (std::vector<int>::iterator c = rand_cells.begin();
         c != rand_cells.end(); c++) {
      // Take into account only those cell neighborhoods for which
      // the central cell contains a catalyzer particle
      cell = local_cells.cell[*c];
      p_local = cell->part;
      np = cell->n;

      // We find all catalyzers in a cell and then randomize their ids
      // for the same reason as above and then start the catalytic
      // reaction procedure
      catalyzers.clear();
      for (int i = 0; i < np; i++) {
        if (p_local[i].p.type == reaction.catalyzer_type)
          catalyzers.push_back(i);
      }
      std::shuffle(catalyzers.begin(), catalyzers.end(), rng);

      // Loop cell neighbors
      for (int n = 0; n < cells.size(); n++) {
        cell = &cells[n];
        p_neigh = cell->part;
        np = cell->n;

        // We loop over all the catalyzer particles
        for (std::vector<int>::iterator id = catalyzers.begin();
             id != catalyzers.end(); id++) {
          reactants.clear();
          products.clear();

          // Particle list loop
          for (int i = 0; i < np; i++) {
            // Get the distance between a catalyst and another particle
            get_mi_vector(vec21, p_local[*id].r.p, p_neigh[i].r.p);
            dist2 = sqrlen(vec21);

            // Check if the distance is within the reaction range and
            // check if no reaction has taken place on the particle in
            // the current step
            if (dist2 < reaction.range * reaction.range &&
                p_neigh[i].p.catalyzer_count == 0) {
              // If the particle is of correct type AND resides in the
              // correct half space, append it to the lists of viable
              // reaction candidates
              if (p_neigh[i].p.type == reaction.reactant_type &&
                  in_lower_half_space(p_local[*id], p_neigh[i])) {
                reactants.push_back(i);
#ifdef ELECTROSTATICS
                reactant_q = p_neigh[i].p.q;
#endif // ELECTROSTATICS
              }
              if (p_neigh[i].p.type == reaction.product_type &&
                  !in_lower_half_space(p_local[*id], p_neigh[i]))
                products.push_back(i);
#ifdef ELECTROSTATICS
              product_q = p_neigh[i].p.q;
#endif // ELECTROSTATICS
            }
          }

          // If reactants and products were found, perform the reaction
          if (reactants.size() > 0 && products.size() > 0) {
            // There cannot be more reactions than the minimum of
            // the number of reactants and products.  Hence we need
            // to determine which number is smaller and also count
            // the number of reactions.
            n_reactions = 0;

            // If there are more products than reactants...
            if (reactants.size() <= products.size()) {
              // ...iterate the reactant...
              for (std::vector<int>::iterator rt = reactants.begin();
                   rt < reactants.end(); rt++) {
                // ...draw a random number number and compare to the
                // reaction rate...
                rand = d_random();
                if (rand > ct_ratexp) {
                  // ...tag the particle for modification...
                  p_neigh[*rt].p.catalyzer_count = 1;
                  n_reactions++;
                }
              }

              // ...tag as many products as there will be reactions
              // at random
              std::shuffle(products.begin(), products.end(), rng);
              for (int p = 0; p < n_reactions; p++)
                p_neigh[products[p]].p.catalyzer_count = 1;
            } else {
              // Same as above, but for the case that the number of
              // reactants is greater than the number of products
              for (std::vector<int>::iterator pt = products.begin();
                   pt < products.end(); pt++) {
                rand = d_random();
                if (rand > ct_ratexp) {
                  p_neigh[*pt].p.catalyzer_count = 1;
                  n_reactions++;
                }
              }

              std::shuffle(reactants.begin(), reactants.end(), rng);
              for (int p = 0; p < n_reactions; p++)
                p_neigh[reactants[p]].p.catalyzer_count = 1;
            }
          }
        }
      }
    }

    // Apply the changes to the tagged particles.  Therefore, again
    // loop over all cells
    for (std::vector<int>::iterator c = rand_cells.begin();
         c != rand_cells.end(); c++) {
      cell = local_cells.cell[*c];
      p_local = cell->part;
      np = cell->n;
      // Particle list loop
      for (int i = 0; i < np; i++) {
        // If the particle has been tagged we perform the changes
        if (p_local[i].p.catalyzer_count != 0) {
          // Flip type and charge
          if (p_local[i].p.type == reaction.reactant_type) {
            p_local[i].p.type = reaction.product_type;
#ifdef ELECTROSTATICS
            p_local[i].p.q = product_q;
#endif // ELECTROSTATICS
          } else {
            p_local[i].p.type = reaction.reactant_type;
#ifdef ELECTROSTATICS
            p_local[i].p.q = reactant_q;
#endif // ELECTROSTATICS
          }
          // Reset the tag for the next step
          p_local[i].p.catalyzer_count = 0;
        }
      }
    }

    on_particle_change();
  }
}
#endif // ROTATION

void integrate_reaction() {
#ifdef ROTATION
  if (reaction.swap)
    integrate_reaction_swap();
  else
#endif // ROTATION
    integrate_reaction_noswap();
}

#endif
