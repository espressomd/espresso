/*
  Copyright (C) 2010-2018 The ESPResSo project
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
/** \file
    Implementation of the energy calculation.
*/
#ifndef ENERGY_INLINE_HPP
#define ENERGY_INLINE_HPP

#include "config.hpp"

#include "bonded_interactions/angle_cosine.hpp"
#include "bonded_interactions/angle_cossquare.hpp"
#include "bonded_interactions/angle_dist.hpp"
#include "bonded_interactions/angle_harmonic.hpp"
#include "bonded_interactions/bonded_interaction_data.hpp"
#include "bonded_interactions/bonded_tab.hpp"
#include "bonded_interactions/dihedral.hpp"
#include "bonded_interactions/fene.hpp"
#include "bonded_interactions/harmonic.hpp"
#include "bonded_interactions/harmonic_dumbbell.hpp"
#include "bonded_interactions/quartic.hpp"
#include "bonded_interactions/subt_lj.hpp"
#include "bonded_interactions/thermalized_bond.hpp"
#include "bonded_interactions/umbrella.hpp"
#include "electrostatics_magnetostatics/debye_hueckel.hpp"
#include "nonbonded_interactions/bmhtf-nacl.hpp"
#include "nonbonded_interactions/buckingham.hpp"
#include "nonbonded_interactions/gaussian.hpp"
#include "nonbonded_interactions/gb.hpp"
#include "nonbonded_interactions/hat.hpp"
#include "nonbonded_interactions/hertzian.hpp"
#include "nonbonded_interactions/lj.hpp"
#include "nonbonded_interactions/ljcos.hpp"
#include "nonbonded_interactions/ljcos2.hpp"
#include "nonbonded_interactions/ljgen.hpp"
#include "nonbonded_interactions/morse.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "nonbonded_interactions/nonbonded_tab.hpp"
#include "nonbonded_interactions/reaction_field.hpp"
#include "nonbonded_interactions/soft_sphere.hpp"
#include "nonbonded_interactions/steppot.hpp"
#include "nonbonded_interactions/thole.hpp"
#ifdef ELECTROSTATICS
#include "bonded_interactions/bonded_coulomb.hpp"
#endif
#ifdef P3M
#include "bonded_interactions/bonded_coulomb_p3m_sr.hpp"
#endif
#include "electrostatics_magnetostatics/elc.hpp"
#include "electrostatics_magnetostatics/mmm1d.hpp"
#include "electrostatics_magnetostatics/mmm2d.hpp"
#include "electrostatics_magnetostatics/p3m-dipolar.hpp"
#include "electrostatics_magnetostatics/p3m.hpp"
#include "electrostatics_magnetostatics/scafacos.hpp"
#include "statistics.hpp"
#include "thermostat.hpp"

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
inline double calc_non_bonded_pair_energy(const Particle *p1,
                                          const Particle *p2,
                                          const IA_parameters *ia_params,
                                          const double d[3], double dist,
                                          double dist2) {
  double ret = 0;

#ifdef NO_INTRA_NB
  if (p1->e->p.mol_id == p2->e->p.mol_id)
    return 0;
#endif

#ifdef LENNARD_JONES
  /* Lennard-Jones */
  ret += lj_pair_energy(p1, p2, ia_params, d, dist);
#endif

#ifdef LENNARD_JONES_GENERIC
  /* Generic Lennard-Jones */
  ret += ljgen_pair_energy(p1, p2, ia_params, d, dist);
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
  /* Morse */
  ret += morse_pair_energy(p1, p2, ia_params, d, dist);
#endif

#ifdef BUCKINGHAM
  /* Buckingham */
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
  /* Lennard-Jones */
  ret += ljcos2_pair_energy(p1, p2, ia_params, d, dist);
#endif

#ifdef THOLE
  /* Thole damping */
  ret += thole_pair_energy(p1, p2, ia_params, d, dist);
#endif

#ifdef TABULATED
  /* tabulated */
  ret += tabulated_pair_energy(p1, p2, ia_params, d, dist);
#endif

#ifdef LJCOS
  /* Lennard-Jones cosine */
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

/** Add non bonded energies and short range Coulomb between a pair of particles.
    @param p1        pointer to particle 1.
    @param p2        pointer to particle 2.
    @param d         vector between p1 and p2.
    @param dist      distance between p1 and p2.
    @param dist2     distance squared between p1 and p2.
*/
inline void add_non_bonded_pair_energy(Particle *p1, Particle *p2, double d[3],
                                       double dist, double dist2) {
  IA_parameters *ia_params = get_ia_param(p1->e->p.type, p2->e->p.type);

#if defined(ELECTROSTATICS) || defined(DIPOLES)
  double ret = 0;
#endif

#ifdef EXCLUSIONS
  if (do_nonbonded(p1, p2))
#endif
    *obsstat_nonbonded(&energy, p1->e->p.type, p2->e->p.type) +=
        calc_non_bonded_pair_energy(p1, p2, ia_params, d, dist, dist2);

#ifdef ELECTROSTATICS
  if (coulomb.method != COULOMB_NONE) {
    /* real space Coulomb */
    switch (coulomb.method) {
#ifdef P3M
    case COULOMB_P3M_GPU:
    case COULOMB_P3M:
      ret = p3m_pair_energy(p1->e->p.q * p2->e->p.q, dist);
      break;
    case COULOMB_ELC_P3M:
      ret = p3m_pair_energy(p1->e->p.q * p2->e->p.q, dist);
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
      ret = mmm2d_coulomb_pair_energy(p1->e->p.q * p2->e->p.q, d, dist2, dist);
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

/** Calculate bonded energies for one particle.
    @param p1 particle for which to calculate energies
*/

inline void add_bonded_energy(Particle *p1) {
  Particle *p3 = nullptr, *p4 = nullptr;
  Bonded_ia_parameters *iaparams;
  int i, bond_broken;
  double ret = 0, dx[3] = {0, 0, 0};

  i = 0;
  while (i < p1->bl.n) {
    int type_num = p1->bl.e[i++];
    iaparams = &bonded_ia_params[type_num];
    int type = iaparams->type;
    int n_partners = iaparams->num;

    /* fetch particle 2, which is always needed */
    Particle *p2 = local_particles[p1->bl.e[i++]];
    if (!p2) {
      runtimeErrorMsg() << "bond broken between particles " << p1->e->p.identity
                        << " and " << p1->bl.e[i - 1]
                        << " (particles not stored on the same node)";
      return;
    }

    /* fetch particle 3 eventually */
    if (n_partners >= 2) {
      p3 = local_particles[p1->bl.e[i++]];
      if (!p3) {
        runtimeErrorMsg() << "bond broken between particles " << p1->e->p.identity
                          << ", " << p1->bl.e[i - 2] << " and "
                          << p1->bl.e[i - 1]
                          << " (particles not stored on the same node)";
        return;
      }
    }

    /* fetch particle 4 eventually */
    if (n_partners >= 3) {
      p4 = local_particles[p1->bl.e[i++]];
      if (!p4) {
        runtimeErrorMsg() << "bond broken between particles " << p1->e->p.identity
                          << ", " << p1->bl.e[i - 3] << ", " << p1->bl.e[i - 2]
                          << " and " << p1->bl.e[i - 1]
                          << " (particles not stored on the same node)";
        return;
      }
    }
    /* similar to the force, we prepare the center-center vector */
    if (n_partners == 1)
      get_mi_vector(dx, p1->r.p, p2->r.p);

    switch (type) {
    case BONDED_IA_FENE:
      bond_broken = fene_pair_energy(p1, p2, iaparams, dx, &ret);
      break;
#ifdef ROTATION
    case BONDED_IA_HARMONIC_DUMBBELL:
      bond_broken = harmonic_dumbbell_pair_energy(p1, p2, iaparams, dx, &ret);
      break;
#endif
    case BONDED_IA_HARMONIC:
      bond_broken = harmonic_pair_energy(p1, p2, iaparams, dx, &ret);
      break;
    case BONDED_IA_QUARTIC:
      bond_broken = quartic_pair_energy(p1, p2, iaparams, dx, &ret);
      break;
    case BONDED_IA_THERMALIZED_DIST:
      bond_broken = thermalized_bond_energy(p1, p2, iaparams, dx, &ret);
      break;
#ifdef ELECTROSTATICS
    case BONDED_IA_BONDED_COULOMB:
      bond_broken = bonded_coulomb_pair_energy(p1, p2, iaparams, dx, &ret);
      break;
#endif
#ifdef P3M
    case BONDED_IA_BONDED_COULOMB_P3M_SR:
      bond_broken =
          bonded_coulomb_p3m_sr_pair_energy(p1, p2, iaparams, dx, &ret);
      break;
#endif
#ifdef LENNARD_JONES
    case BONDED_IA_SUBT_LJ:
      bond_broken = subt_lj_pair_energy(p1, p2, iaparams, dx, &ret);
      break;
#endif
#ifdef BOND_ANGLE
    case BONDED_IA_ANGLE_HARMONIC:
      bond_broken = angle_harmonic_energy(p1, p2, p3, iaparams, &ret);
      break;
    case BONDED_IA_ANGLE_COSINE:
      bond_broken = angle_cosine_energy(p1, p2, p3, iaparams, &ret);
      break;
    case BONDED_IA_ANGLE_COSSQUARE:
      bond_broken = angle_cossquare_energy(p1, p2, p3, iaparams, &ret);
      break;
#endif
#ifdef BOND_ANGLEDIST
    case BONDED_IA_ANGLEDIST:
      bond_broken = angledist_energy(p1, p2, p3, iaparams, &ret);
      break;
#endif
    case BONDED_IA_DIHEDRAL:
      bond_broken = dihedral_energy(p2, p1, p3, p4, iaparams, &ret);
      break;
#ifdef BOND_CONSTRAINT
    case BONDED_IA_RIGID_BOND:
      bond_broken = 0;
      ret = 0;
      break;
#endif
#ifdef TABULATED
    case BONDED_IA_TABULATED:
      switch (iaparams->p.tab.type) {
      case TAB_BOND_LENGTH:
        bond_broken = tab_bond_energy(p1, p2, iaparams, dx, &ret);
        break;
      case TAB_BOND_ANGLE:
        bond_broken = tab_angle_energy(p1, p2, p3, iaparams, &ret);
        break;
      case TAB_BOND_DIHEDRAL:
        bond_broken = tab_dihedral_energy(p2, p1, p3, p4, iaparams, &ret);
        break;
      default:
        runtimeErrorMsg() << "add_bonded_energy: tabulated bond type of atom "
                          << p1->e->p.identity << " unknown\n";
        return;
      }
      break;
#endif
#ifdef UMBRELLA
    case BONDED_IA_UMBRELLA:
      bond_broken = umbrella_pair_energy(p1, p2, iaparams, dx, &ret);
      break;
#endif
    case BONDED_IA_VIRTUAL_BOND:
      bond_broken = 0;
      ret = 0;
      break;
    default:
      runtimeErrorMsg() << "add_bonded_energy: bond type (" << type
                        << ") of atom " << p1->e->p.identity << " unknown\n";
      return;
    }

    if (bond_broken) {
      switch (n_partners) {
      case 1: {
        runtimeErrorMsg() << "bond broken between particles " << p1->e->p.identity
                          << " and " << p2->e->p.identity;
        break;
      }
      case 2: {
        runtimeErrorMsg() << "bond broken between particles " << p1->e->p.identity
                          << ", " << p2->e->p.identity << " and "
                          << p3->e->p.identity;
        break;
      }
      case 3: {
        runtimeErrorMsg() << "bond broken between particles " << p1->e->p.identity
                          << ", " << p2->e->p.identity << ", " << p3->e->p.identity
                          << " and " << p4->e->p.identity;
        break;
      }
      }
      // bond broken, don't add whatever we find in the energy
      continue;
    }

    *obsstat_bonded(&energy, type_num) += ret;
  }
}

/** Calculate kinetic energies for one particle.
    @param p1 particle for which to calculate energies
*/
inline void add_kinetic_energy(Particle *p1) {
#ifdef VIRTUAL_SITES
  if (p1->e->p.is_virtual)
    return;
#endif

  /* kinetic energy */
  energy.data.e[0] += (Utils::sqr(p1->e->m.v[0]) + Utils::sqr(p1->e->m.v[1]) +
                       Utils::sqr(p1->e->m.v[2])) *
                      0.5 * p1->e->p.mass;

#ifdef ROTATION
  if (p1->e->p.rotation) {
    /* the rotational part is added to the total kinetic energy;
       Here we use the rotational inertia  */

    energy.data.e[0] += 0.5 * (Utils::sqr(p1->e->m.omega[0]) * p1->e->p.rinertia[0] +
                               Utils::sqr(p1->e->m.omega[1]) * p1->e->p.rinertia[1] +
                               Utils::sqr(p1->e->m.omega[2]) * p1->e->p.rinertia[2]);
  }
#endif
}

inline void add_single_particle_energy(Particle *p) {
  add_kinetic_energy(p);
  add_bonded_energy(p);
}

#endif // ENERGY_INLINE_HPP
