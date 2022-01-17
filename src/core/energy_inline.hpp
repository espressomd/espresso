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

#include "energy.hpp"

#include "bonded_interactions/bonded_interaction_data.hpp"
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
#include "electrostatics_magnetostatics/coulomb_inline.hpp"
#endif

#ifdef DIPOLES
#include "electrostatics_magnetostatics/dipole_inline.hpp"
#endif

#include "Observable_stat.hpp"
#include "Particle.hpp"
#include "bond_error.hpp"
#include "errorhandling.hpp"
#include "exclusions.hpp"

#include <utils/Span.hpp>
#include <utils/Vector.hpp>

#include <boost/optional.hpp>
#include <boost/range/algorithm/find_if.hpp>
#include <boost/variant.hpp>

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
 *  to the energy observable.
 *  @param p1        particle 1.
 *  @param p2        particle 2.
 *  @param d         vector between p1 and p2.
 *  @param dist      distance between p1 and p2.
 *  @param dist2     distance squared between p1 and p2.
 *  @param[in,out] obs_energy   energy observable.
 */
inline void add_non_bonded_pair_energy(Particle const &p1, Particle const &p2,
                                       Utils::Vector3d const &d,
                                       double const dist, double const dist2,
                                       Observable_stat &obs_energy) {
  IA_parameters const &ia_params = *get_ia_param(p1.p.type, p2.p.type);

#ifdef EXCLUSIONS
  if (do_nonbonded(p1, p2))
#endif
    obs_energy.add_non_bonded_contribution(
        p1.p.type, p2.p.type,
        calc_non_bonded_pair_energy(p1, p2, ia_params, d, dist));

#ifdef ELECTROSTATICS
  if (!obs_energy.coulomb.empty())
    obs_energy.coulomb[0] +=
        Coulomb::pair_energy(p1, p2, p1.p.q * p2.p.q, d, dist, dist2);
#endif

#ifdef DIPOLES
  if (!obs_energy.dipolar.empty())
    obs_energy.dipolar[0] += Dipole::pair_energy(p1, p2, d, dist, dist2);
#endif
}

inline boost::optional<double>
calc_bonded_energy(Bonded_IA_Parameters const &iaparams, Particle const &p1,
                   Utils::Span<Particle *> partners) {
  auto const n_partners = partners.size();

  auto p2 = (n_partners > 0) ? partners[0] : nullptr;
  auto p3 = (n_partners > 1) ? partners[1] : nullptr;
  auto p4 = (n_partners > 2) ? partners[2] : nullptr;

  if (n_partners == 1) {
    auto const dx = box_geo.get_mi_vector(p1.r.p, p2->r.p);
    if (auto const *iap = boost::get<FeneBond>(&iaparams)) {
      return iap->energy(dx);
    }
    if (auto const *iap = boost::get<HarmonicBond>(&iaparams)) {
      return iap->energy(dx);
    }
    if (auto const *iap = boost::get<QuarticBond>(&iaparams)) {
      return iap->energy(dx);
    }
#ifdef ELECTROSTATICS
    if (auto const *iap = boost::get<BondedCoulomb>(&iaparams)) {
      return iap->energy(p1.p.q * p2->p.q, dx);
    }
    if (auto const *iap = boost::get<BondedCoulombSR>(&iaparams)) {
      return iap->energy(p1, *p2, dx);
    }
#endif
#ifdef BOND_CONSTRAINT
    if (auto const *iap = boost::get<RigidBond>(&iaparams)) {
      return boost::optional<double>(0);
    }
#endif
#ifdef TABULATED
    if (auto const *iap = boost::get<TabulatedDistanceBond>(&iaparams)) {
      return iap->energy(dx);
    }
#endif
    if (auto const *iap = boost::get<VirtualBond>(&iaparams)) {
      return boost::optional<double>(0);
    }
    throw BondUnknownTypeError();
  } // 1 partner
  if (n_partners == 2) {
    if (auto const *iap = boost::get<AngleHarmonicBond>(&iaparams)) {
      return iap->energy(p1.r.p, p2->r.p, p3->r.p);
    }
    if (auto const *iap = boost::get<AngleCosineBond>(&iaparams)) {
      return iap->energy(p1.r.p, p2->r.p, p3->r.p);
    }
    if (auto const *iap = boost::get<AngleCossquareBond>(&iaparams)) {
      return iap->energy(p1.r.p, p2->r.p, p3->r.p);
    }
    if (auto const *iap = boost::get<TabulatedAngleBond>(&iaparams)) {
      return iap->energy(p1.r.p, p2->r.p, p3->r.p);
    }
    if (auto const *iap = boost::get<IBMTriel>(&iaparams)) {
      runtimeWarningMsg() << "Unsupported bond type " +
                                 std::to_string(iaparams.which()) +
                                 " in energy calculation.";
      return 0.;
    }
    throw BondUnknownTypeError();
  } // 2 partners
  if (n_partners == 3) {
    if (auto const *iap = boost::get<DihedralBond>(&iaparams)) {
      return iap->energy(p2->r.p, p1.r.p, p3->r.p, p4->r.p);
    }
    if (auto const *iap = boost::get<TabulatedDihedralBond>(&iaparams)) {
      return iap->energy(p2->r.p, p1.r.p, p3->r.p, p4->r.p);
    }
    if (auto const *iap = boost::get<IBMTribend>(&iaparams)) {
      runtimeWarningMsg() << "Unsupported bond type " +
                                 std::to_string(iaparams.which()) +
                                 " in energy calculation.";
      return 0.;
    }
    throw BondUnknownTypeError();
  } // 3 partners
  if (n_partners == 0) {
    return 0.;
  }

  throw BondInvalidSizeError(n_partners);
}

/** Calculate kinetic energies from translation for one particle.
 *  @param p   particle for which to calculate energies
 */
inline double translational_kinetic_energy(Particle const &p) {
  return p.p.is_virtual ? 0. : 0.5 * p.p.mass * p.m.v.norm2();
}

/** Calculate kinetic energies from rotation for one particle.
 *  @param p   particle for which to calculate energies
 */
inline double rotational_kinetic_energy(Particle const &p) {
#ifdef ROTATION
  return p.p.rotation
             ? 0.5 * (hadamard_product(p.m.omega, p.m.omega) * p.p.rinertia)
             : 0.0;
#else
  return 0.0;
#endif
}

/** Calculate kinetic energies for one particle.
 *  @param p   particle for which to calculate energies
 */
inline double calc_kinetic_energy(Particle const &p) {
  return translational_kinetic_energy(p) + rotational_kinetic_energy(p);
}

#endif // ENERGY_INLINE_HPP
