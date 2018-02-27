/*
  Copyright (C) 2010,2012,2013,2014,2015,2016,2017 The ESPResSo project
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
/** \file pressure_inline.hpp
    Pressure calculation. Really similar to \ref energy.hpp "energy.h".
*/

#ifndef CORE_PRESSURE_INLINE_HPP
#define CORE_PRESSURE_INLINE_HPP

#include "debug.hpp"
#include "forces_inline.hpp"
#include "integrate.hpp"
#include "npt.hpp"
#include "pressure.hpp"
#include "thermostat.hpp"
#include "utils.hpp"

/** Calculate non bonded energies between a pair of particles.
    @param p1        pointer to particle 1.
    @param p2        pointer to particle 2.
    @param d         vector between p1 and p2.
    @param dist      distance between p1 and p2.
    @param dist2     distance squared between p1 and p2. */
inline void add_non_bonded_pair_virials(Particle *p1, Particle *p2, double d[3],
                                        double dist, double dist2) {
  int p1molid, p2molid, k, l;
  double force[3] = {0, 0, 0};

#ifdef EXCLUSIONS
  if (do_nonbonded(p1, p2))
#endif
  {
    calc_non_bonded_pair_force(p1, p2, d, dist, dist2, force);
    *obsstat_nonbonded(&virials, p1->p.type, p2->p.type) +=
        d[0] * force[0] + d[1] * force[1] + d[2] * force[2];

    /* stress tensor part */
    for (k = 0; k < 3; k++)
      for (l = 0; l < 3; l++)
        obsstat_nonbonded(&p_tensor, p1->p.type, p2->p.type)[k * 3 + l] +=
            force[k] * d[l];

    p1molid = p1->p.mol_id;
    p2molid = p2->p.mol_id;
    if (p1molid == p2molid) {
      *obsstat_nonbonded_intra(&virials_non_bonded, p1->p.type, p2->p.type) +=
          d[0] * force[0] + d[1] * force[1] + d[2] * force[2];

      for (k = 0; k < 3; k++)
        for (l = 0; l < 3; l++)
          obsstat_nonbonded_intra(&p_tensor_non_bonded, p1->p.type,
                                  p2->p.type)[k * 3 + l] += force[k] * d[l];
    }
    if (p1molid != p2molid) {
      *obsstat_nonbonded_inter(&virials_non_bonded, p1->p.type, p2->p.type) +=
          d[0] * force[0] + d[1] * force[1] + d[2] * force[2];

      for (k = 0; k < 3; k++)
        for (l = 0; l < 3; l++)
          obsstat_nonbonded_inter(&p_tensor_non_bonded, p1->p.type,
                                  p2->p.type)[k * 3 + l] += force[k] * d[l];
    }
  }

#ifdef ELECTROSTATICS
  /* real space coulomb */
  if (coulomb.method != COULOMB_NONE) {
    switch (coulomb.method) {
#ifdef P3M
    case COULOMB_P3M_GPU:
    case COULOMB_P3M:
      force[0] = 0.0;
      force[1] = 0.0;
      force[2] = 0.0;
      p3m_add_pair_force(p1->p.q * p2->p.q, d, dist2, dist, force);
      virials.coulomb[0] += p3m_pair_energy(p1->p.q * p2->p.q, dist);
      for (k = 0; k < 3; k++)
        for (l = 0; l < 3; l++)
          p_tensor.coulomb[k * 3 + l] += force[k] * d[l];

      break;
#endif

    /* short range potentials, where we use the virial */
    /***************************************************/
    case COULOMB_DH: {
      double force[3] = {0, 0, 0};

      add_dh_coulomb_pair_force(p1, p2, d, dist, force);
      for (k = 0; k < 3; k++)
        for (l = 0; l < 3; l++)
          p_tensor.coulomb[k * 3 + l] += force[k] * d[l];
      virials.coulomb[0] += force[0] * d[0] + force[1] * d[1] + force[2] * d[2];
      break;
    }
    case COULOMB_RF: {
      double force[3] = {0, 0, 0};

      add_rf_coulomb_pair_force(p1, p2, d, dist, force);
      for (k = 0; k < 3; k++)
        for (l = 0; l < 3; l++)
          p_tensor.coulomb[k * 3 + l] += force[k] * d[l];
      virials.coulomb[0] += force[0] * d[0] + force[1] * d[1] + force[2] * d[2];
      break;
    }
    case COULOMB_INTER_RF:
      // this is done together with the other short range interactions
      break;
    default:
      fprintf(stderr, "calculating pressure for electrostatics method that "
                      "doesn't have it implemented\n");
      break;
    }
  }
#endif /*ifdef ELECTROSTATICS */

#ifdef DIPOLES
  /* real space magnetic dipole-dipole */
  if (coulomb.Dmethod != DIPOLAR_NONE) {
    fprintf(stderr, "calculating pressure for magnetostatics which doesn't "
                    "have it implemented\n");
  }
#endif /*ifdef DIPOLES */
}

/*New add_bonded_virials function*/
// only pair bonds contribute in this function
#ifdef BOND_CLASS_DEBUG
inline void add_bonded_virials(Particle *p1)
{

  int bond_broken = bond_container.virial_loop(p1);
  if(bond_broken == 2){
    return;
  };

}
#endif // BOND_CLASS_DEBUG


/*New add_bonded_virials function*/
// only pair bonds contribute in this function
#ifdef BOND_CLASS_DEBUG
inline void add_three_body_bonded_stress(Particle *p1)
{

  int bond_broken = bond_container.three_body_stress_loop(p1);
  if(bond_broken == 2){
    return;
  };

}
#endif // BOND_CLASS_DEBUG


/** Calculate kinetic pressure (aka energy) for one particle.
    @param p1 particle for which to calculate pressure
    @param v_comp flag which enables (1) compensation of the velocities required
                  for deriving a pressure reflecting \ref nptiso_struct::p_inst
                  (hence it only works with domain decomposition); naturally it
                  therefore doesn't make sense to use it without NpT.
*/
inline void add_kinetic_virials(Particle *p1, int v_comp) {
  int k, l;
  /* kinetic energy */
#ifdef MULTI_TIMESTEP
  if (smaller_time_step > 0.) {
    if (v_comp)
      virials.data.e[0] += Utils::sqr(time_step / smaller_time_step) *
                           (Utils::sqr(p1->m.v[0] - p1->f.f[0]) +
                            Utils::sqr(p1->m.v[1] - p1->f.f[1]) +
                            Utils::sqr(p1->m.v[2] - p1->f.f[2])) *
                           (*p1).p.mass;
    else
      virials.data.e[0] += Utils::sqr(time_step / smaller_time_step) *
                           (Utils::sqr(p1->m.v[0]) + Utils::sqr(p1->m.v[1]) +
                            Utils::sqr(p1->m.v[2])) *
                           (*p1).p.mass;
  } else
#endif
  {
    if (v_comp)
      virials.data.e[0] += (Utils::sqr(p1->m.v[0] - p1->f.f[0]) +
                            Utils::sqr(p1->m.v[1] - p1->f.f[1]) +
                            Utils::sqr(p1->m.v[2] - p1->f.f[2])) *
                           (*p1).p.mass;
    else
      virials.data.e[0] += (Utils::sqr(p1->m.v[0]) + Utils::sqr(p1->m.v[1]) +
                            Utils::sqr(p1->m.v[2])) *
                           (*p1).p.mass;
  }

  /* ideal gas contribution (the rescaling of the velocities by '/=time_step'
   * each will be done later) */
  for (k = 0; k < 3; k++)
    for (l = 0; l < 3; l++)
#ifdef MULTI_TIMESTEP
      if (smaller_time_step > 0.)
        p_tensor.data.e[k * 3 + l] +=
            Utils::sqr(time_step / smaller_time_step) * (p1->m.v[k]) *
            (p1->m.v[l]) * (*p1).p.mass;
      else
#endif
        p_tensor.data.e[k * 3 + l] +=
            (p1->m.v[k]) * (p1->m.v[l]) * (*p1).p.mass;
}

#endif
