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
 *  Implementation of the energy calculation.
 */
#ifndef ENERGY_INLINE_HPP
#define ENERGY_INLINE_HPP

#include "config.hpp"

#include "bonded_interactions/angle_cosine.hpp"
#include "bonded_interactions/angle_cossquare.hpp"
#include "bonded_interactions/angle_harmonic.hpp"
#include "bonded_interactions/bonded_interaction_data.hpp"
#include "bonded_interactions/bonded_tab.hpp"
#include "bonded_interactions/dihedral.hpp"
#include "bonded_interactions/fene.hpp"
#include "bonded_interactions/harmonic.hpp"
#include "bonded_interactions/harmonic_dumbbell.hpp"
#include "bonded_interactions/quartic.hpp"
#include "bonded_interactions/subt_lj.hpp"
#include "bonded_interactions/umbrella.hpp"
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
#include "nonbonded_interactions/soft_sphere.hpp"
#include "nonbonded_interactions/steppot.hpp"
#include "nonbonded_interactions/thole.hpp"
#include "nonbonded_interactions/wca.hpp"
#ifdef ELECTROSTATICS
#include "bonded_interactions/bonded_coulomb.hpp"
#include "bonded_interactions/bonded_coulomb_sr.hpp"
#include "electrostatics_magnetostatics/coulomb_inline.hpp"
#endif
#include "statistics.hpp"

#include "energy.hpp"

#ifdef DIPOLES
#include "electrostatics_magnetostatics/dipole_inline.hpp"
#endif

/** Calculate non bonded energies between a pair of particles.
 *  @param p1         pointer to particle 1.
 *  @param p2         pointer to particle 2.
 *  @param ia_params  the interaction parameters between the two particles
 *  @param d          vector between p1 and p2.
 *  @param dist       distance between p1 and p2.
 *  @param dist2      distance squared between p1 and p2.
 *  @return the short ranged interaction energy between the two particles
 */
inline double calc_non_bonded_pair_energy(const Particle *p1,
                                          const Particle *p2,
                                          const IA_parameters *ia_params,
                                          const double d[3], double dist,
                                          double dist2) {
  double ret = 0;

#ifdef NO_INTRA_NB
  if (p1->p.mol_id == p2->p.mol_id)
    return 0;
#endif

#ifdef LENNARD_JONES
  /* Lennard-Jones */
  ret += lj_pair_energy(ia_params, dist);
#endif
#ifdef WCA
  /* WCA */
  ret += wca_pair_energy(p1, p2, ia_params, d, dist);
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

  return ret;
}

/** Add non bonded energies and short range Coulomb between a pair of particles.
 *  @param p1        pointer to particle 1.
 *  @param p2        pointer to particle 2.
 *  @param d         vector between p1 and p2.
 *  @param dist      distance between p1 and p2.
 *  @param dist2     distance squared between p1 and p2.
 */
inline void add_non_bonded_pair_energy(const Particle *p1, const Particle *p2,
                                       const double *d, double dist,
                                       double dist2) {
  IA_parameters const *ia_params = get_ia_param(p1->p.type, p2->p.type);

#if defined(ELECTROSTATICS) || defined(DIPOLES)
#endif

#ifdef EXCLUSIONS
  if (do_nonbonded(p1, p2))
#endif
    *obsstat_nonbonded(&energy, p1->p.type, p2->p.type) +=
        calc_non_bonded_pair_energy(p1, p2, ia_params, d, dist, dist2);

#ifdef ELECTROSTATICS
  energy.coulomb[0] +=
      Coulomb::pair_energy(p1, p2, p1->p.q * p2->p.q, d, dist, dist2);
#endif

#ifdef DIPOLES
  Dipole::add_pair_energy(p1, p2, d, dist, dist2, energy);
#endif
}

/** Calculate bonded energies for one particle.
 *  @param p1 particle for which to calculate energies
 */
inline void add_bonded_energy(const Particle *p1) {
  Particle *p3 = nullptr, *p4 = nullptr;
  Bonded_ia_parameters *iaparams;
  int i, bond_broken = 1;
  double ret = 0;

  i = 0;
  while (i < p1->bl.n) {
    int type_num = p1->bl.e[i++];
    iaparams = &bonded_ia_params[type_num];
    int type = iaparams->type;
    int n_partners = iaparams->num;

    /* fetch particle 2, which is always needed */
    Particle *p2 = local_particles[p1->bl.e[i++]];
    if (!p2) {
      runtimeErrorMsg() << "bond broken between particles " << p1->p.identity
                        << " and " << p1->bl.e[i - 1]
                        << " (particles not stored on the same node)";
      return;
    }

    /* fetch particle 3 eventually */
    if (n_partners >= 2) {
      p3 = local_particles[p1->bl.e[i++]];
      if (!p3) {
        runtimeErrorMsg() << "bond broken between particles " << p1->p.identity
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
        runtimeErrorMsg() << "bond broken between particles " << p1->p.identity
                          << ", " << p1->bl.e[i - 3] << ", " << p1->bl.e[i - 2]
                          << " and " << p1->bl.e[i - 1]
                          << " (particles not stored on the same node)";
        return;
      }
    }

    if (n_partners == 1) {
      auto const dx = get_mi_vector(p1->r.p, p2->r.p);
      switch (type) {
      case BONDED_IA_FENE:
        bond_broken = fene_pair_energy(iaparams, dx, &ret);
        break;
#ifdef ROTATION
      case BONDED_IA_HARMONIC_DUMBBELL:
        bond_broken = harmonic_dumbbell_pair_energy(p1, iaparams, dx, &ret);
        break;
#endif
      case BONDED_IA_HARMONIC:
        bond_broken = harmonic_pair_energy(iaparams, dx, &ret);
        break;
      case BONDED_IA_QUARTIC:
        bond_broken = quartic_pair_energy(iaparams, dx, &ret);
        break;
#ifdef ELECTROSTATICS
      case BONDED_IA_BONDED_COULOMB:
        bond_broken = bonded_coulomb_pair_energy(p1, p2, iaparams, dx, &ret);
        break;
      case BONDED_IA_BONDED_COULOMB_SR:
        bond_broken = bonded_coulomb_sr_pair_energy(p1, p2, iaparams, dx, &ret);
        break;
#endif
#ifdef LENNARD_JONES
      case BONDED_IA_SUBT_LJ:
        bond_broken = subt_lj_pair_energy(p1, p2, iaparams, dx, &ret);
        break;
#endif
#ifdef BOND_CONSTRAINT
      case BONDED_IA_RIGID_BOND:
        bond_broken = 0;
        ret = 0;
        break;
#endif
#ifdef TABULATED
      case BONDED_IA_TABULATED:
        if (iaparams->num == 1)
          bond_broken = tab_bond_energy(iaparams, dx, &ret);
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
                          << ") of atom " << p1->p.identity << " unknown\n";
        return;
      }
    } // 1 partner
    else if (n_partners == 2) {
      switch (type) {
      case BONDED_IA_ANGLE_HARMONIC:
        bond_broken = angle_harmonic_energy(p1, p2, p3, iaparams, &ret);
        break;
      case BONDED_IA_ANGLE_COSINE:
        bond_broken = angle_cosine_energy(p1, p2, p3, iaparams, &ret);
        break;
      case BONDED_IA_ANGLE_COSSQUARE:
        bond_broken = angle_cossquare_energy(p1, p2, p3, iaparams, &ret);
        break;
#ifdef TABULATED
      case BONDED_IA_TABULATED:
        if (iaparams->num == 2)
          bond_broken = tab_angle_energy(p1, p2, p3, iaparams, &ret);
        break;
#endif
      default:
        runtimeErrorMsg() << "add_bonded_energy: bond type (" << type
                          << ") of atom " << p1->p.identity << " unknown\n";
        return;
      }
    } // 2 partner
    else if (n_partners == 3) {
      switch (type) {
      case BONDED_IA_DIHEDRAL:
        bond_broken = dihedral_energy(p2, p1, p3, p4, iaparams, &ret);
        break;
#ifdef TABULATED
      case BONDED_IA_TABULATED:
        if (iaparams->num == 3)
          bond_broken = tab_dihedral_energy(p1, p2, p3, p4, iaparams, &ret);
        break;
#endif
      default:
        runtimeErrorMsg() << "add_bonded_energy: bond type (" << type
                          << ") of atom " << p1->p.identity << " unknown\n";
        return;
      }
    } // 3 partners

    if (bond_broken) {
      switch (n_partners) {
      case 1: {
        runtimeErrorMsg() << "bond broken between particles " << p1->p.identity
                          << " and " << p2->p.identity;
        break;
      }
      case 2: {
        runtimeErrorMsg() << "bond broken between particles " << p1->p.identity
                          << ", " << p2->p.identity << " and "
                          << p3->p.identity;
        break;
      }
      case 3: {
        runtimeErrorMsg() << "bond broken between particles " << p1->p.identity
                          << ", " << p2->p.identity << ", " << p3->p.identity
                          << " and " << p4->p.identity;
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
 *  @param p1 particle for which to calculate energies
 */
inline void add_kinetic_energy(const Particle *p1) {
#ifdef VIRTUAL_SITES
  if (p1->p.is_virtual)
    return;
#endif

  /* kinetic energy */
  energy.data.e[0] += 0.5 * p1->p.mass * p1->m.v.norm2();

#ifdef ROTATION
  if (p1->p.rotation) {
    /* the rotational part is added to the total kinetic energy;
       Here we use the rotational inertia  */

    energy.data.e[0] += 0.5 * (Utils::sqr(p1->m.omega[0]) * p1->p.rinertia[0] +
                               Utils::sqr(p1->m.omega[1]) * p1->p.rinertia[1] +
                               Utils::sqr(p1->m.omega[2]) * p1->p.rinertia[2]);
  }
#endif
}

inline void add_single_particle_energy(const Particle *p) {
  add_kinetic_energy(p);
  add_bonded_energy(p);
}

#endif // ENERGY_INLINE_HPP
