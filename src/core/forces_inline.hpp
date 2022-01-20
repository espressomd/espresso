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
#ifndef _FORCES_INLINE_HPP
#define _FORCES_INLINE_HPP

#include "config.hpp"

#include "forces.hpp"

#include "bond_breakage.hpp"
#include "bonded_interactions/bonded_interaction_data.hpp"
#include "bonded_interactions/thermalized_bond_kernel.hpp"
#include "immersed_boundary/ibm_tribend.hpp"
#include "immersed_boundary/ibm_triel.hpp"
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
#include "object-in-fluid/oif_global_forces.hpp"
#include "object-in-fluid/oif_local_forces.hpp"

#ifdef DIPOLES
#include "electrostatics_magnetostatics/dipole_inline.hpp"
#endif

#ifdef ELECTROSTATICS
#include "electrostatics_magnetostatics/coulomb_inline.hpp"
#endif

#ifdef DPD
#include "dpd.hpp"
#endif

#include "Particle.hpp"
#include "bond_error.hpp"
#include "errorhandling.hpp"
#include "exclusions.hpp"
#include "thermostat.hpp"

#include <utils/Span.hpp>
#include <utils/Vector.hpp>

#include <boost/optional.hpp>
#include <boost/variant.hpp>

#include <tuple>

inline ParticleForce calc_non_bonded_pair_force(Particle const &p1,
                                                Particle const &p2,
                                                IA_parameters const &ia_params,
                                                Utils::Vector3d const &d,
                                                double const dist) {

  ParticleForce pf{};
  double force_factor = 0;
/* Lennard-Jones */
#ifdef LENNARD_JONES
  force_factor += lj_pair_force_factor(ia_params, dist);
#endif
/* WCA */
#ifdef WCA
  force_factor += wca_pair_force_factor(ia_params, dist);
#endif
/* Lennard-Jones generic */
#ifdef LENNARD_JONES_GENERIC
  force_factor += ljgen_pair_force_factor(ia_params, dist);
#endif
/* smooth step */
#ifdef SMOOTH_STEP
  force_factor += SmSt_pair_force_factor(ia_params, dist);
#endif
/* Hertzian force */
#ifdef HERTZIAN
  force_factor += hertzian_pair_force_factor(ia_params, dist);
#endif
/* Gaussian force */
#ifdef GAUSSIAN
  force_factor += gaussian_pair_force_factor(ia_params, dist);
#endif
/* BMHTF NaCl */
#ifdef BMHTF_NACL
  force_factor += BMHTF_pair_force_factor(ia_params, dist);
#endif
/* Buckingham*/
#ifdef BUCKINGHAM
  force_factor += buck_pair_force_factor(ia_params, dist);
#endif
/* Morse*/
#ifdef MORSE
  force_factor += morse_pair_force_factor(ia_params, dist);
#endif
/*soft-sphere potential*/
#ifdef SOFT_SPHERE
  force_factor += soft_pair_force_factor(ia_params, dist);
#endif
/*hat potential*/
#ifdef HAT
  force_factor += hat_pair_force_factor(ia_params, dist);
#endif
/* Lennard-Jones cosine */
#ifdef LJCOS
  force_factor += ljcos_pair_force_factor(ia_params, dist);
#endif
/* Lennard-Jones cosine */
#ifdef LJCOS2
  force_factor += ljcos2_pair_force_factor(ia_params, dist);
#endif
/* Thole damping */
#ifdef THOLE
  pf.f += thole_pair_force(p1, p2, ia_params, d, dist);
#endif
/* tabulated */
#ifdef TABULATED
  force_factor += tabulated_pair_force_factor(ia_params, dist);
#endif
/* Gay-Berne */
#ifdef GAY_BERNE
  // The gb force function isn't inlined, probably due to its size
  if (dist < ia_params.gay_berne.cut) {
    pf += gb_pair_force(p1.r.calc_director(), p2.r.calc_director(), ia_params,
                        d, dist);
  }
#endif
  pf.f += force_factor * d;
  return pf;
}

inline ParticleForce calc_opposing_force(ParticleForce const &pf,
                                         Utils::Vector3d const &d) {
  ParticleForce out{-pf.f};
#ifdef ROTATION
  // if torque is a null vector, the opposing torque is a null vector too
  // (this check guards from returning a small yet non-null opposing
  // torque due to numerical imprecision)
  if (pf.torque[0] != 0. || pf.torque[1] != 0. || pf.torque[2] != 0.) {
    out.torque = -(pf.torque + vector_product(d, pf.f));
  }
#endif
  return out;
}

/** Calculate non-bonded forces between a pair of particles and update their
 *  forces and torques.
 *  @param[in,out] p1   particle 1.
 *  @param[in,out] p2   particle 2.
 *  @param[in] d        vector between @p p1 and @p p2.
 *  @param dist         distance between @p p1 and @p p2.
 *  @param dist2        distance squared between @p p1 and @p p2.
 */
inline void add_non_bonded_pair_force(Particle &p1, Particle &p2,
                                      Utils::Vector3d const &d, double dist,
                                      double dist2) {
  IA_parameters const &ia_params = *get_ia_param(p1.p.type, p2.p.type);
  ParticleForce pf{};

  /***********************************************/
  /* non-bonded pair potentials                  */
  /***********************************************/

  if (dist < ia_params.max_cut) {
#ifdef EXCLUSIONS
    if (do_nonbonded(p1, p2))
#endif
      pf += calc_non_bonded_pair_force(p1, p2, ia_params, d, dist);
  }

  /***********************************************/
  /* short-range electrostatics                  */
  /***********************************************/

#ifdef ELECTROSTATICS
  {
    auto const forces = Coulomb::pair_force(p1, p2, d, dist);
    pf.f += std::get<0>(forces);
#ifdef P3M
    // forces from the virtual charges
    p1.f.f += std::get<1>(forces);
    p2.f.f += std::get<2>(forces);
#endif
  }
#endif

  /*********************************************************************/
  /* everything before this contributes to the virial pressure in NpT, */
  /* but nothing afterwards                                            */
  /*********************************************************************/
#ifdef NPT
  npt_add_virial_force_contribution(pf.f, d);
#endif

  /***********************************************/
  /* thermostat                                  */
  /***********************************************/

  /* The inter dpd force should not be part of the virial */
#ifdef DPD
  if (thermo_switch & THERMO_DPD) {
    auto const force = dpd_pair_force(p1, p2, ia_params, d, dist, dist2);
    p1.f.f += force;
    p2.f.f -= force;
  }
#endif

  /***********************************************/
  /* long-range magnetostatics                   */
  /***********************************************/

#ifdef DIPOLES
  /* real space magnetic dipole-dipole */
  pf += Dipole::pair_force(p1, p2, d, dist, dist2);
#endif

  /***********************************************/
  /* add total non-bonded forces to particles    */
  /***********************************************/

  p1.f += pf;
  p2.f += calc_opposing_force(pf, d);
}

/** Compute the bonded interaction force between particle pairs.
 *
 *  @param[in] p1          First particle.
 *  @param[in] p2          Second particle.
 *  @param[in] iaparams    Bonded parameters for the interaction.
 *  @param[in] dx          Vector between @p p1 and @p p2.
 */
inline boost::optional<Utils::Vector3d>
calc_bond_pair_force(Particle const &p1, Particle const &p2,
                     Bonded_IA_Parameters const &iaparams,
                     Utils::Vector3d const &dx) {
  if (auto const *iap = boost::get<FeneBond>(&iaparams)) {
    return iap->force(dx);
  }
  if (auto const *iap = boost::get<HarmonicBond>(&iaparams)) {
    return iap->force(dx);
  }
  if (auto const *iap = boost::get<QuarticBond>(&iaparams)) {
    return iap->force(dx);
  }
#ifdef ELECTROSTATICS
  if (auto const *iap = boost::get<BondedCoulomb>(&iaparams)) {
    return iap->force(p1.p.q * p2.p.q, dx);
  }
  if (auto const *iap = boost::get<BondedCoulombSR>(&iaparams)) {
    return iap->force(dx);
  }
#endif
#ifdef TABULATED
  if (auto const *iap = boost::get<TabulatedDistanceBond>(&iaparams)) {
    return iap->force(dx);
  }
#endif
  if (boost::get<VirtualBond>(&iaparams) || boost::get<RigidBond>(&iaparams)) {
    return Utils::Vector3d{};
  }
  throw BondUnknownTypeError();
}

inline bool add_bonded_two_body_force(Bonded_IA_Parameters const &iaparams,
                                      Particle &p1, Particle &p2) {
  auto const dx = box_geo.get_mi_vector(p1.r.p, p2.r.p);

  if (auto const *iap = boost::get<ThermalizedBond>(&iaparams)) {
    auto result = iap->forces(p1, p2, dx);
    if (result) {
      using std::get;
      p1.f.f += get<0>(result.get());
      p2.f.f += get<1>(result.get());

      return false;
    }
  } else {
    auto result = calc_bond_pair_force(p1, p2, iaparams, dx);
    if (result) {
      p1.f.f += result.get();
      p2.f.f -= result.get();

#ifdef NPT
      npt_add_virial_force_contribution(result.get(), dx);
#endif
      return false;
    }
  }
  return true;
}

inline boost::optional<
    std::tuple<Utils::Vector3d, Utils::Vector3d, Utils::Vector3d>>
calc_bonded_three_body_force(Bonded_IA_Parameters const &iaparams,
                             Particle const &p1, Particle const &p2,
                             Particle const &p3) {
  if (auto const *iap = boost::get<AngleHarmonicBond>(&iaparams)) {
    return iap->forces(p1.r.p, p2.r.p, p3.r.p);
  }
  if (auto const *iap = boost::get<AngleCosineBond>(&iaparams)) {
    return iap->forces(p1.r.p, p2.r.p, p3.r.p);
  }
  if (auto const *iap = boost::get<AngleCossquareBond>(&iaparams)) {
    return iap->forces(p1.r.p, p2.r.p, p3.r.p);
  }
#ifdef TABULATED
  if (auto const *iap = boost::get<TabulatedAngleBond>(&iaparams)) {
    return iap->forces(p1.r.p, p2.r.p, p3.r.p);
  }
#endif
  if (auto const *iap = boost::get<IBMTriel>(&iaparams)) {
    return iap->calc_forces(p1, p2, p3);
  }
  throw BondUnknownTypeError();
}

inline bool add_bonded_three_body_force(Bonded_IA_Parameters const &iaparams,
                                        Particle &p1, Particle &p2,
                                        Particle &p3) {
  if (auto const *iap = boost::get<OifGlobalForcesBond>(&iaparams)) {
    return false;
  }
  auto const result = calc_bonded_three_body_force(iaparams, p1, p2, p3);
  if (result) {
    using std::get;
    auto const &forces = result.get();

    p1.f.f += get<0>(forces);
    p2.f.f += get<1>(forces);
    p3.f.f += get<2>(forces);

    return false;
  }
  return true;
}

inline boost::optional<std::tuple<Utils::Vector3d, Utils::Vector3d,
                                  Utils::Vector3d, Utils::Vector3d>>
calc_bonded_four_body_force(Bonded_IA_Parameters const &iaparams,
                            Particle const &p1, Particle const &p2,
                            Particle const &p3, Particle const &p4) {
  if (auto const *iap = boost::get<OifLocalForcesBond>(&iaparams)) {
    return iap->calc_forces(p1, p2, p3, p4);
  }
  if (auto const *iap = boost::get<IBMTribend>(&iaparams)) {
    return iap->calc_forces(p1, p2, p3, p4);
  }
  if (auto const *iap = boost::get<DihedralBond>(&iaparams)) {
    return iap->forces(p2.r.p, p1.r.p, p3.r.p, p4.r.p);
  }
#ifdef TABULATED
  if (auto const *iap = boost::get<TabulatedDihedralBond>(&iaparams)) {
    return iap->forces(p2.r.p, p1.r.p, p3.r.p, p4.r.p);
  }
#endif
  throw BondUnknownTypeError();
}

inline bool add_bonded_four_body_force(Bonded_IA_Parameters const &iaparams,
                                       Particle &p1, Particle &p2, Particle &p3,
                                       Particle &p4) {
  auto const result = calc_bonded_four_body_force(iaparams, p1, p2, p3, p4);
  if (result) {
    using std::get;
    auto const &forces = result.get();

    p1.f.f += get<0>(forces);
    p2.f.f += get<1>(forces);
    p3.f.f += get<2>(forces);
    p4.f.f += get<3>(forces);

    return false;
  }

  return true;
}

inline bool add_bonded_force(Particle &p1, int bond_id,
                             Utils::Span<Particle *> partners) {

  // Consider for bond breakage
  if (partners.size() == 1) {
    auto d = box_geo.get_mi_vector(p1.r.p, partners[0]->r.p).norm();
    if (BondBreakage::check_and_handle_breakage(
            p1.p.identity, partners[0]->p.identity, bond_id, d))
      return false;
  }

  auto const &iaparams = *bonded_ia_params.at(bond_id);

  switch (number_of_partners(iaparams)) {
  case 0:
    return false;
  case 1:
    return add_bonded_two_body_force(iaparams, p1, *partners[0]);
  case 2:
    return add_bonded_three_body_force(iaparams, p1, *partners[0],
                                       *partners[1]);
  case 3:
    return add_bonded_four_body_force(iaparams, p1, *partners[0], *partners[1],
                                      *partners[2]);
  default:
    throw BondInvalidSizeError{number_of_partners(iaparams)};
  }
}

#endif
