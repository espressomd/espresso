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
#include "errorhandling.hpp"
#include "nonbonded_interactions/bmhtf-nacl.hpp"
#include "nonbonded_interactions/buckingham.hpp"
#include "nonbonded_interactions/gaussian.hpp"
#include "nonbonded_interactions/gay_berne.hpp"
#include "nonbonded_interactions/hat.hpp"
#include "nonbonded_interactions/hertzian.hpp"
#include "nonbonded_interactions/lj.hpp"
#include "nonbonded_interactions/ljcos.hpp"
#include "nonbonded_interactions/ljcos2.hpp"
#include "nonbonded_interactions/ljgen.hpp"
#include "nonbonded_interactions/morse.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "nonbonded_interactions/nonbonded_tab.hpp"
#include "nonbonded_interactions/smooth_step.hpp"
#include "nonbonded_interactions/soft_sphere.hpp"
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

/** Calculate non-bonded energies between a pair of particles.
 *  @param p1         particle 1.
 *  @param p2         particle 2.
 *  @param ia_params  the interaction parameters between the two particles
 *  @param d          vector between p1 and p2.
 *  @param dist       distance between p1 and p2.
 *  @return the short-range interaction energy between the two particles
 */
inline double calc_non_bonded_pair_energy(Particle const &p1,
                                          Particle const &p2,
                                          IA_parameters const &ia_params,
                                          Utils::Vector3d const &d,
                                          double const dist) {
#ifdef NO_INTRA_NB
  if (p1.p.mol_id == p2.p.mol_id)
    return 0;
#endif

  double ret = 0;

#ifdef LENNARD_JONES
  /* Lennard-Jones */
  ret += lj_pair_energy(ia_params, dist);
#endif
#ifdef WCA
  /* WCA */
  ret += wca_pair_energy(ia_params, dist);
#endif

#ifdef LENNARD_JONES_GENERIC
  /* Generic Lennard-Jones */
  ret += ljgen_pair_energy(ia_params, dist);
#endif

#ifdef SMOOTH_STEP
  /* smooth step */
  ret += SmSt_pair_energy(ia_params, dist);
#endif

#ifdef HERTZIAN
  /* Hertzian potential */
  ret += hertzian_pair_energy(ia_params, dist);
#endif

#ifdef GAUSSIAN
  /* Gaussian potential */
  ret += gaussian_pair_energy(ia_params, dist);
#endif

#ifdef BMHTF_NACL
  /* BMHTF NaCl */
  ret += BMHTF_pair_energy(ia_params, dist);
#endif

#ifdef MORSE
  /* Morse */
  ret += morse_pair_energy(ia_params, dist);
#endif

#ifdef BUCKINGHAM
  /* Buckingham */
  ret += buck_pair_energy(ia_params, dist);
#endif

#ifdef SOFT_SPHERE
  /* soft-sphere */
  ret += soft_pair_energy(ia_params, dist);
#endif

#ifdef HAT
  /* hat */
  ret += hat_pair_energy(ia_params, dist);
#endif

#ifdef LJCOS2
  /* Lennard-Jones */
  ret += ljcos2_pair_energy(ia_params, dist);
#endif

#ifdef THOLE
  /* Thole damping */
  ret += thole_pair_energy(p1, p2, ia_params, d, dist);
#endif

#ifdef TABULATED
  /* tabulated */
  ret += tabulated_pair_energy(ia_params, dist);
#endif

#ifdef LJCOS
  /* Lennard-Jones cosine */
  ret += ljcos_pair_energy(ia_params, dist);
#endif

#ifdef GAY_BERNE
  /* Gay-Berne */
  ret += gb_pair_energy(p1.r.calc_director(), p2.r.calc_director(), ia_params,
                        d, dist);
#endif

  return ret;
}

/** Add non-bonded and short-range Coulomb energies between a pair of particles
 *  to the @ref energy observable.
 *  @param p1        particle 1.
 *  @param p2        particle 2.
 *  @param d         vector between p1 and p2.
 *  @param dist      distance between p1 and p2.
 *  @param dist2     distance squared between p1 and p2.
 */
inline void add_non_bonded_pair_energy(Particle const &p1, Particle const &p2,
                                       Utils::Vector3d const &d,
                                       double const dist, double const dist2) {
  IA_parameters const &ia_params = *get_ia_param(p1.p.type, p2.p.type);

#ifdef EXCLUSIONS
  if (do_nonbonded(p1, p2))
#endif
    *obsstat_nonbonded(&energy, p1.p.type, p2.p.type) +=
        calc_non_bonded_pair_energy(p1, p2, ia_params, d, dist);

#ifdef ELECTROSTATICS
  energy.coulomb[0] +=
      Coulomb::pair_energy(p1, p2, p1.p.q * p2.p.q, d, dist, dist2);
#endif

#ifdef DIPOLES
  energy.dipolar[0] += Dipole::pair_energy(p1, p2, d, dist, dist2);
#endif
}

/** Add bonded energies for one particle to the @ref energy observable.
 *  @param[in] p1   particle for which to calculate energies
 */
inline void add_bonded_energy(Particle const *const p1) {
  int i = 0;
  while (i < p1->bl.n) {
    Particle const *p3 = nullptr;
    Particle const *p4 = nullptr;
    int type_num = p1->bl.e[i++];
    Bonded_ia_parameters const &iaparams = bonded_ia_params[type_num];
    int type = iaparams.type;
    int n_partners = iaparams.num;

    /* fetch particle 2, which is always needed */
    Particle const *const p2 = local_particles[p1->bl.e[i++]];
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

    boost::optional<double> retval;

    if (n_partners == 1) {
      auto const dx = get_mi_vector(p1->r.p, p2->r.p, box_geo);
      switch (type) {
      case BONDED_IA_FENE:
        retval = fene_pair_energy(iaparams, dx);
        break;
#ifdef ROTATION
      case BONDED_IA_HARMONIC_DUMBBELL:
        retval =
            harmonic_dumbbell_pair_energy(p1->r.calc_director(), iaparams, dx);
        break;
#endif
      case BONDED_IA_HARMONIC:
        retval = harmonic_pair_energy(iaparams, dx);
        break;
      case BONDED_IA_QUARTIC:
        retval = quartic_pair_energy(iaparams, dx);
        break;
#ifdef ELECTROSTATICS
      case BONDED_IA_BONDED_COULOMB:
        retval = bonded_coulomb_pair_energy(p1->p.q * p2->p.q, iaparams, dx);
        break;
      case BONDED_IA_BONDED_COULOMB_SR:
        retval = bonded_coulomb_sr_pair_energy(*p1, *p2, iaparams, dx);
        break;
#endif
#ifdef LENNARD_JONES
      case BONDED_IA_SUBT_LJ:
        retval = subt_lj_pair_energy(*get_ia_param(p1->p.type, p2->p.type), dx);
        break;
#endif
#ifdef BOND_CONSTRAINT
      case BONDED_IA_RIGID_BOND:
        retval = boost::optional<double>(0);
        break;
#endif
      case BONDED_IA_TABULATED_DISTANCE:
        retval = tab_bond_energy(iaparams, dx);
        break;
#ifdef UMBRELLA
      case BONDED_IA_UMBRELLA:
        retval = umbrella_pair_energy(iaparams, dx);
        break;
#endif
      case BONDED_IA_VIRTUAL_BOND:
        retval = boost::optional<double>(0);
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
        retval = boost::optional<double>(
            angle_harmonic_energy(p1->r.p, p2->r.p, p3->r.p, iaparams));
        break;
      case BONDED_IA_ANGLE_COSINE:
        retval = boost::optional<double>(
            angle_cosine_energy(p1->r.p, p2->r.p, p3->r.p, iaparams));
        break;
      case BONDED_IA_ANGLE_COSSQUARE:
        retval = boost::optional<double>(
            angle_cossquare_energy(p1->r.p, p2->r.p, p3->r.p, iaparams));
        break;
      case BONDED_IA_TABULATED_ANGLE:
        retval = boost::optional<double>(
            tab_angle_energy(p1->r.p, p2->r.p, p3->r.p, iaparams));
        break;
      default:
        runtimeErrorMsg() << "add_bonded_energy: bond type (" << type
                          << ") of atom " << p1->p.identity << " unknown\n";
        return;
      }
    } // 2 partner
    else if (n_partners == 3) {
      switch (type) {
      case BONDED_IA_DIHEDRAL:
        retval = dihedral_energy(p2->r.p, p1->r.p, p3->r.p, p4->r.p, iaparams);
        break;
      case BONDED_IA_TABULATED_DIHEDRAL:
        retval =
            tab_dihedral_energy(p2->r.p, p1->r.p, p3->r.p, p4->r.p, iaparams);
        break;
      default:
        runtimeErrorMsg() << "add_bonded_energy: bond type (" << type
                          << ") of atom " << p1->p.identity << " unknown\n";
        return;
      }
    } // 3 partners

    if (retval) {
      *obsstat_bonded(&energy, type_num) += retval.get();
    } else {
      switch (n_partners) {
      case 1:
        runtimeErrorMsg() << "bond broken between particles " << p1->p.identity
                          << " and " << p2->p.identity;
        break;
      case 2:
        runtimeErrorMsg() << "bond broken between particles " << p1->p.identity
                          << ", " << p2->p.identity << " and "
                          << p3->p.identity;
        break;
      case 3:
        runtimeErrorMsg() << "bond broken between particles " << p1->p.identity
                          << ", " << p2->p.identity << ", " << p3->p.identity
                          << " and " << p4->p.identity;
        break;
      }
      // bond broken, don't add whatever we find in the energy
      continue;
    }
  }
}

/** Add kinetic energies for one particle to the @ref energy observable.
 *  @param[in] p1   particle for which to calculate energies
 */
inline void add_kinetic_energy(Particle const &p1) {
#ifdef VIRTUAL_SITES
  if (p1.p.is_virtual)
    return;
#endif

  /* kinetic energy */
  if (not p1.p.is_virtual)
    energy.data.e[0] += 0.5 * p1.p.mass * p1.m.v.norm2();

    // Note that rotational degrees of virtual sites are integrated
    // and therefore can contribute to kinetic energy
#ifdef ROTATION
  if (p1.p.rotation) {
    /* the rotational part is added to the total kinetic energy;
       Here we use the rotational inertia  */

    energy.data.e[0] += 0.5 * (Utils::sqr(p1.m.omega[0]) * p1.p.rinertia[0] +
                               Utils::sqr(p1.m.omega[1]) * p1.p.rinertia[1] +
                               Utils::sqr(p1.m.omega[2]) * p1.p.rinertia[2]);
  }
#endif
}

/** Add kinetic and bonded energies for one particle to the @ref energy
 *  observable.
 *  @param[in] p   particle for which to calculate energies
 */
inline void add_single_particle_energy(Particle const &p) {
  add_kinetic_energy(p);
  add_bonded_energy(&p);
}

#endif // ENERGY_INLINE_HPP
