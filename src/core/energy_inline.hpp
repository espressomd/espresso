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
/** \file energy_inline.hpp
    Implementation of the energy calculation.
*/
#ifndef ENERGY_INLINE_HPP
#define ENERGY_INLINE_HPP

#include "config.hpp"

#include "bmhtf-nacl.hpp"
#include "buckingham.hpp"
#include "dihedral.hpp"
#include "thermalized_bond.hpp"
#include "fene.hpp"
#include "gaussian.hpp"
#include "gb.hpp"
#include "harmonic.hpp"
#include "harmonic_dumbbell.hpp"
#include "hat.hpp"
#include "hertzian.hpp"
#include "lj.hpp"
#include "ljangle.hpp"
#include "ljcos.hpp"
#include "ljcos2.hpp"
#include "ljgen.hpp"
#include "overlap.hpp"
#include "p3m-dipolar.hpp"
#include "p3m.hpp"
#include "quartic.hpp"
#include "soft_sphere.hpp"
#include "statistics.hpp"
#include "steppot.hpp"
#include "tab.hpp"
#include "thole.hpp"
#include "thermostat.hpp"
#include "umbrella.hpp"
#ifdef ELECTROSTATICS
#include "bonded_coulomb.hpp"
#endif
#ifdef P3M
#include "bonded_coulomb_p3m_sr.hpp"
#endif
#include "angle_cosine.hpp"
#include "angle_cossquare.hpp"
#include "angle_harmonic.hpp"
#include "angledist.hpp"
#include "debye_hueckel.hpp"
#include "elc.hpp"
#include "endangledist.hpp"
#include "hydrogen_bond.hpp"
#include "mmm1d.hpp"
#include "mmm2d.hpp"
#include "morse.hpp"
#include "reaction_field.hpp"
#include "scafacos.hpp"
#include "subt_lj.hpp"
#include "twist_stack.hpp"

#ifdef CONSTRAINTS
#include "constraints.hpp"
#endif

#ifdef EXTERNAL_FORCES
#include "external_potential.hpp"
#endif

#include "energy.hpp"

/** Calculate non bonded energies between a pair of particles.
    @param p1        pointer to particle 1.
    @param p2        pointer to particle 2.
    @param ia_params the interaction parameters between the two particles
    @param d         vector between p1 and p2.
    @param dist      distance between p1 and p2.
    @param dist2     distance squared between p1 and p2.
    @return the short ranged interaction energy between the two particles
*/
inline double calc_non_bonded_pair_energy(const Particle *p1, const Particle *p2,
                                          const IA_parameters *ia_params, const double d[3],
                                          double dist, double dist2) {
  double ret = 0;

#ifdef NO_INTRA_NB
  if (p1->p.mol_id == p2->p.mol_id)
    return 0;
#endif

#ifdef LENNARD_JONES
  /* lennard jones */
  ret += lj_pair_energy(p1, p2, ia_params, d, dist);
#endif

#ifdef LENNARD_JONES_GENERIC
  /* Generic lennard jones */
  ret += ljgen_pair_energy(p1, p2, ia_params, d, dist);
#endif

#ifdef LJ_ANGLE
  /* Directional LJ */
  ret += ljangle_pair_energy(p1, p2, ia_params, d, dist);
#endif

#ifdef SMOOTH_STEP
  /* smooth step */
  ret += SmSt_pair_energy(p1, p2, ia_params, d, dist, dist2);
#endif

#ifdef HERTZIAN
  /* Hertzian potential */
  ret += hertzian_pair_energy(p1, p2, ia_params, d, dist, dist2);
#endif

#ifdef GAUSSIAN
  /* Gaussian potential */
  ret += gaussian_pair_energy(p1, p2, ia_params, d, dist, dist2);
#endif

#ifdef BMHTF_NACL
  /* BMHTF NaCl */
  ret += BMHTF_pair_energy(p1, p2, ia_params, d, dist, dist2);
#endif

#ifdef MORSE
  /* morse */
  ret += morse_pair_energy(p1, p2, ia_params, d, dist);
#endif

#ifdef BUCKINGHAM
  /* lennard jones */
  ret += buck_pair_energy(p1, p2, ia_params, d, dist);
#endif

#ifdef SOFT_SPHERE
  /* soft-sphere */
  ret += soft_pair_energy(p1, p2, ia_params, d, dist);
#endif

#ifdef HAT
  /* hat */
  ret += hat_pair_energy(p1, p2, ia_params, d, dist);
#endif

#ifdef LJCOS2
  /* lennard jones */
  ret += ljcos2_pair_energy(p1, p2, ia_params, d, dist);
#endif

#ifdef THOLE
  /* thole damping */
  ret += thole_pair_energy(p1, p2, ia_params, d, dist);
#endif

#ifdef TABULATED
  /* tabulated */
  ret += tabulated_pair_energy(p1, p2, ia_params, d, dist);
#endif

#ifdef LJCOS
  /* lennard jones cosine */
  ret += ljcos_pair_energy(p1, p2, ia_params, d, dist);
#endif

#ifdef GAY_BERNE
  /* Gay-Berne */
  ret += gb_pair_energy(p1, p2, ia_params, d, dist);
#endif

#ifdef INTER_RF
  ret += interrf_pair_energy(p1, p2, ia_params, dist);
#endif

  return ret;
}

/** Add non bonded energies and short range coulomb between a pair of particles.
    @param p1        pointer to particle 1.
    @param p2        pointer to particle 2.
    @param d         vector between p1 and p2.
    @param dist      distance between p1 and p2.
    @param dist2     distance squared between p1 and p2.
*/
inline void add_non_bonded_pair_energy(Particle *p1, Particle *p2, double d[3],
                                       double dist, double dist2) {
  IA_parameters *ia_params = get_ia_param(p1->p.type, p2->p.type);

#if defined(ELECTROSTATICS) || defined(DIPOLES)
  double ret = 0;
#endif

#ifdef EXCLUSIONS
  if (do_nonbonded(p1, p2))
#endif
  *obsstat_nonbonded(&energy, p1->p.type, p2->p.type) +=
      calc_non_bonded_pair_energy(p1, p2, ia_params, d, dist, dist2);

#ifdef ELECTROSTATICS
  if (coulomb.method != COULOMB_NONE) {
    /* real space coulomb */
    switch (coulomb.method) {
#ifdef P3M
    case COULOMB_P3M_GPU:
    case COULOMB_P3M:
      ret = p3m_pair_energy(p1->p.q * p2->p.q, dist);
      break;
    case COULOMB_ELC_P3M:
      ret = p3m_pair_energy(p1->p.q * p2->p.q, dist);
      if (elc_params.dielectric_contrast_on)
        ret += 0.5 * ELC_P3M_dielectric_layers_energy_contribution(p1, p2);
      break;
#endif
#ifdef SCAFACOS
    case COULOMB_SCAFACOS:
      ret += Scafacos::pair_energy(p1, p2, dist);
      break;
#endif
    case COULOMB_DH:
      ret = dh_coulomb_pair_energy(p1, p2, dist);
      break;
    case COULOMB_RF:
      ret = rf_coulomb_pair_energy(p1, p2, dist);
      break;
    case COULOMB_INTER_RF:
      // this is done above as interaction
      ret = 0;
      break;
    case COULOMB_MMM1D:
      ret = mmm1d_coulomb_pair_energy(p1, p2, d, dist2, dist);
      break;
    case COULOMB_MMM2D:
      ret = mmm2d_coulomb_pair_energy(p1->p.q * p2->p.q, d, dist2, dist);
      break;
    default:
      ret = 0.;
    }
    energy.coulomb[0] += ret;
  }
#endif

#ifdef DIPOLES
  if (coulomb.Dmethod != DIPOLAR_NONE) {
    // ret=0;
    switch (coulomb.Dmethod) {
#ifdef DP3M
    case DIPOLAR_MDLC_P3M:
    // fall trough
    case DIPOLAR_P3M:
      ret = dp3m_pair_energy(p1, p2, d, dist2, dist);
      break;
#endif
    default:
      ret = 0;
    }
    energy.dipolar[0] += ret;
  }
#endif
}

#ifdef BOND_CLASS_DEBUG
/**New add_bonded_energy function for bond classes**/
inline void add_bonded_energy(Particle *p1)
{
  int bond_broken = bond_container.energy_loop(p1);
  if(bond_broken == 2)
    return;
}
#endif //BOND_CLASS_DEBUG

/** Calculate kinetic energies for one particle.
    @param p1 particle for which to calculate energies
*/
inline void add_kinetic_energy(Particle *p1) {
#ifdef VIRTUAL_SITES
  if (p1->p.isVirtual)
    return;
#endif

  /* kinetic energy */

  // #ifdef MULTI_TIMESTEP
  //   if (p1->p.smaller_timestep==1) {
  //     ostringstream msg;
  //     msg << "SMALL TIME STEP";
  //     energy.data.e[0] += Utils::sqr(smaller_time_step/time_step) *
  //       (Utils::sqr(p1->m.v[0]) + Utils::sqr(p1->m.v[1]) + Utils::sqr(p1->m.v[2]))*(*p1).p.mass;
  //   }
  //   else
  // #endif
  energy.data.e[0] +=
      (Utils::sqr(p1->m.v[0]) + Utils::sqr(p1->m.v[1]) + Utils::sqr(p1->m.v[2])) * (*p1).p.mass;

#ifdef ROTATION
  if (p1->p.rotation)
  {
    /* the rotational part is added to the total kinetic energy;
       Here we use the rotational inertia  */

    energy.data.e[0] += (Utils::sqr(p1->m.omega[0]) * p1->p.rinertia[0] +
                         Utils::sqr(p1->m.omega[1]) * p1->p.rinertia[1] +
                         Utils::sqr(p1->m.omega[2]) * p1->p.rinertia[2]) *
                        time_step * time_step;
  }
#endif
}

inline void add_single_particle_energy(Particle *p) {
  add_kinetic_energy(p);
  add_bonded_energy(p);
#ifdef CONSTRAINTS
  add_constraints_energy(p);
#endif
#ifdef EXTERNAL_FORCES
  add_external_potential_energy(p);
#endif
}

#endif // ENERGY_INLINE_HPP
