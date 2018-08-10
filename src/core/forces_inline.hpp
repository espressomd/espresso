/*
  Copyright (C) 2010,2012,2013,2014,2015,2016 The ESPResSo project
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
#ifndef _FORCES_INLINE_HPP
#define _FORCES_INLINE_HPP

#include "config.hpp"

#ifdef MOLFORCES
#include "topology.hpp"
#endif

#include "thermostat.hpp"
#include "electrostatics/mmm1d.hpp"
#include "electrostatics/mmm2d.hpp"
#include "electrostatics/p3m.hpp"
#include "bonded_interactions/angle_cosine.hpp"
#include "bonded_interactions/angle_cossquare.hpp"
#include "bonded_interactions/angle_harmonic.hpp"
#include "bonded_interactions/angle_dist.hpp"
#include "bonded_interactions/quartic.hpp"
#include "bonded_interactions/dihedral.hpp"
#include "bonded_interactions/thermalized_bond.hpp"
#include "bonded_interactions/fene.hpp"
#include "bonded_interactions/harmonic.hpp"
#include "bonded_interactions/harmonic_dumbbell.hpp"
#include "bonded_interactions/subt_lj.hpp"
#include "bonded_interactions/umbrella.hpp"
#include "bonded_interactions/bonded_tab.hpp"
#include "nonbonded_interactions/nonbonded_tab.hpp"
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
#include "nonbonded_interactions/thole.hpp"
#include "nonbonded_interactions/morse.hpp"
#include "nonbonded_interactions/soft_sphere.hpp"
#include "nonbonded_interactions/steppot.hpp"
#include "magnetic_non_p3m_methods.hpp"
#include "mdlc_correction.hpp"
#include "metadynamics.hpp"
#include "electrostatics/mmm1d.hpp"
#include "electrostatics/mmm2d.hpp"
#include "molforces.hpp"
#include "npt.hpp"
#include "object-in-fluid/affinity.hpp"
#include "object-in-fluid/membrane_collision.hpp"
#include "object-in-fluid/oif_global_forces.hpp"
#include "object-in-fluid/oif_local_forces.hpp"
#include "object-in-fluid/out_direction.hpp"
#include "electrostatics/p3m-dipolar.hpp"
#include "electrostatics/p3m.hpp"
#include "collision.hpp"
#include "electrostatics/elc.hpp"
#include "forces.hpp"
#ifdef ELECTROSTATICS
#include "bonded_interactions/bonded_coulomb.hpp"
#include "electrostatics/debye_hueckel.hpp"
#include "nonbonded_interactions/reaction_field.hpp"
#include "electrostatics/scafacos.hpp"
#endif
#ifdef P3M
#include "bonded_interactions/bonded_coulomb_p3m_sr.hpp"
#endif
#ifdef IMMERSED_BOUNDARY
#include "immersed_boundary/ibm_triel.hpp"
#include "immersed_boundary/ibm_tribend.hpp"
#endif
#ifdef DPD
#include "dpd.hpp"
#endif

/** initialize the forces for a ghost particle */
inline void init_ghost_force(Particle *part) {
  part->f.f[0] = 0;
  part->f.f[1] = 0;
  part->f.f[2] = 0;

#ifdef ROTATION
  {
    double scale;
    /* set torque to zero */
    part->f.torque[0] = 0;
    part->f.torque[1] = 0;
    part->f.torque[2] = 0;

    /* and rescale quaternion, so it is exactly of unit length */
    scale = sqrt(Utils::sqr(part->r.quat[0]) + Utils::sqr(part->r.quat[1]) +
                 Utils::sqr(part->r.quat[2]) + Utils::sqr(part->r.quat[3]));
    if (scale == 0) {
      part->r.quat[0] = 1;
    } else {
      part->r.quat[0] /= scale;
      part->r.quat[1] /= scale;
      part->r.quat[2] /= scale;
      part->r.quat[3] /= scale;
    }
  }
#endif
}

/** initialize the forces for a real particle */
inline void init_local_particle_force(Particle *part) {
  if (thermo_switch & THERMO_LANGEVIN)
    friction_thermo_langevin(part);
  else {
    part->f.f[0] = 0;
    part->f.f[1] = 0;
    part->f.f[2] = 0;
  }

#ifdef EXTERNAL_FORCES
  if (part->p.ext_flag & PARTICLE_EXT_FORCE) {
    part->f.f[0] += part->p.ext_force[0];
    part->f.f[1] += part->p.ext_force[1];
    part->f.f[2] += part->p.ext_force[2];
  }
#endif

#ifdef ROTATION
  {
    double scale;
    /* set torque to zero */
    part->f.torque[0] = 0;
    part->f.torque[1] = 0;
    part->f.torque[2] = 0;

#ifdef EXTERNAL_FORCES
    if (part->p.ext_flag & PARTICLE_EXT_TORQUE) {
      part->f.torque[0] += part->p.ext_torque[0];
      part->f.torque[1] += part->p.ext_torque[1];
      part->f.torque[2] += part->p.ext_torque[2];
    }
#endif

#ifdef ENGINE
    // apply a swimming force in the direction of
    // the particle's orientation axis
    if (part->swim.swimming) {
      part->f.f[0] += part->swim.f_swim * part->r.quatu[0];
      part->f.f[1] += part->swim.f_swim * part->r.quatu[1];
      part->f.f[2] += part->swim.f_swim * part->r.quatu[2];
    }
#endif

    /* and rescale quaternion, so it is exactly of unit length */
    scale = sqrt(Utils::sqr(part->r.quat[0]) + Utils::sqr(part->r.quat[1]) +
                 Utils::sqr(part->r.quat[2]) + Utils::sqr(part->r.quat[3]));
    if (scale == 0) {
      part->r.quat[0] = 1;
    } else {
      part->r.quat[0] /= scale;
      part->r.quat[1] /= scale;
      part->r.quat[2] /= scale;
      part->r.quat[3] /= scale;
    }
  }
#endif
}

inline void calc_non_bonded_pair_force_parts(
    const Particle *const p1, const Particle *const p2,
    IA_parameters *ia_params, double d[3], double dist, double dist2,
    double force[3], double torque1[3] = nullptr, double torque2[3] = nullptr) {
#ifdef NO_INTRA_NB
  if (p1->p.mol_id == p2->p.mol_id)
    return;
#endif
/* lennard jones */
#ifdef LENNARD_JONES
  add_lj_pair_force(p1, p2, ia_params, d, dist, force);
#endif
/* lennard jones generic */
#ifdef LENNARD_JONES_GENERIC
  add_ljgen_pair_force(p1, p2, ia_params, d, dist, force);
#endif
/* smooth step */
#ifdef SMOOTH_STEP
  add_SmSt_pair_force(p1, p2, ia_params, d, dist, dist2, force);
#endif
/* Hertzian force */
#ifdef HERTZIAN
  add_hertzian_pair_force(p1, p2, ia_params, d, dist, dist2, force);
#endif
/* Gaussian force */
#ifdef GAUSSIAN
  add_gaussian_pair_force(p1, p2, ia_params, d, dist, dist2, force);
#endif
/* BMHTF NaCl */
#ifdef BMHTF_NACL
  add_BMHTF_pair_force(p1, p2, ia_params, d, dist, dist2, force);
#endif
/* buckingham*/
#ifdef BUCKINGHAM
  add_buck_pair_force(p1, p2, ia_params, d, dist, force);
#endif
/* morse*/
#ifdef MORSE
  add_morse_pair_force(p1, p2, ia_params, d, dist, force);
#endif
/*soft-sphere potential*/
#ifdef SOFT_SPHERE
  add_soft_pair_force(p1, p2, ia_params, d, dist, force);
#endif
/*repulsive membrane potential*/
#ifdef MEMBRANE_COLLISION
  add_membrane_collision_pair_force(p1, p2, ia_params, d, dist, force);
#endif
/*hat potential*/
#ifdef HAT
  add_hat_pair_force(p1, p2, ia_params, d, dist, force);
#endif
/* lennard jones cosine */
#ifdef LJCOS
  add_ljcos_pair_force(p1, p2, ia_params, d, dist, force);
#endif
/* lennard jones cosine */
#ifdef LJCOS2
  add_ljcos2_pair_force(p1, p2, ia_params, d, dist, force);
#endif
/* thole damping */
#ifdef THOLE
  add_thole_pair_force(p1, p2, ia_params, d, dist, force);
#endif
/* tabulated */
#ifdef TABULATED
  add_tabulated_pair_force(p1, p2, ia_params, d, dist, force);
#endif
/* Gay-Berne */
#ifdef GAY_BERNE
  add_gb_pair_force(p1, p2, ia_params, d, dist, force, torque1, torque2);
#endif
#ifdef INTER_RF
  add_interrf_pair_force(p1, p2, ia_params, d, dist, force);
#endif
}

inline void calc_non_bonded_pair_force(const Particle *p1, const Particle *p2,
                                       IA_parameters *ia_params, double *d,
                                       double dist, double dist2, double *force,
                                       double *torque1 = nullptr,
                                       double *torque2 = nullptr) {
  calc_non_bonded_pair_force_parts(p1, p2, ia_params, d, dist, dist2, force,
                                   torque1, torque2);
}

inline void calc_non_bonded_pair_force(Particle *p1, Particle *p2, double d[3],
                                       double dist, double dist2,
                                       double force[3]) {
  IA_parameters *ia_params = get_ia_param(p1->p.type, p2->p.type);
  calc_non_bonded_pair_force(p1, p2, ia_params, d, dist, dist2, force);
}

/** Calculate non bonded forces between a pair of particles.
    @param p1        pointer to particle 1.
    @param p2        pointer to particle 2.
    @param d         vector between p1 and p2.
    @param dist      distance between p1 and p2.
    @param dist2     distance squared between p1 and p2. */
inline void add_non_bonded_pair_force(Particle *p1, Particle *p2, double d[3],
                                      double dist, double dist2) {

  IA_parameters *ia_params = get_ia_param(p1->p.type, p2->p.type);
  Vector3d force{};
  double torque1[3] = {0., 0., 0.};
  double torque2[3] = {0., 0., 0.};
  int j;

  /***********************************************/
  /* bond creation and breaking                  */
  /***********************************************/

#ifdef COLLISION_DETECTION
  if (collision_params.mode != COLLISION_MODE_OFF)
    detect_collision(p1, p2, dist);
#endif

/*affinity potential*/
#ifdef AFFINITY
  add_affinity_pair_force(p1, p2, ia_params, d, dist, force.data());
#endif

  FORCE_TRACE(fprintf(stderr, "%d: interaction %d<->%d dist %f\n", this_node,
                      p1->p.identity, p2->p.identity, dist));

  /***********************************************/
  /* non bonded pair potentials                  */
  /***********************************************/

#ifdef EXCLUSIONS
  if (do_nonbonded(p1, p2))
#endif
    calc_non_bonded_pair_force(p1, p2, ia_params, d, dist, dist2, force.data(),
                               torque1, torque2);

    /***********************************************/
    /* short range electrostatics                  */
    /***********************************************/

#ifdef ELECTROSTATICS
  if (coulomb.method == COULOMB_DH)
    add_dh_coulomb_pair_force(p1, p2, d, dist, force.data());

  if (coulomb.method == COULOMB_RF)
    add_rf_coulomb_pair_force(p1, p2, d, dist, force.data());
#endif

/*********************************************************************/
/* everything before this contributes to the virial pressure in NpT, */
/* but nothing afterwards                                            */
/*********************************************************************/
#ifdef NPT
  for (j = 0; j < 3; j++)
    if (integ_switch == INTEG_METHOD_NPT_ISO)
      nptiso.p_vir[j] += force[j] * d[j];
#endif

/***********************************************/
/* thermostat                                  */
/***********************************************/

/** The inter dpd force should not be part of the virial */
#ifdef DPD
  if (thermo_switch & THERMO_DPD) {
    force += dpd_pair_force(p1, p2, ia_params, d, dist, dist2);
  }
#endif

  /***********************************************/
  /* long range electrostatics                   */
  /***********************************************/

#ifdef ELECTROSTATICS
  /* real space coulomb */
  const double q1q2 = p1->p.q * p2->p.q;

  switch (coulomb.method) {
#ifdef P3M
  case COULOMB_ELC_P3M: {
    if (q1q2) {
      p3m_add_pair_force(q1q2, d, dist2, dist, force.data());

      // forces from the virtual charges
      // they go directly onto the particles, since they are not pairwise forces
      if (elc_params.dielectric_contrast_on)
        ELC_P3M_dielectric_layers_force_contribution(p1, p2, p1->f.f.data(),
                                                     p2->f.f.data());
    }
    break;
  }
  case COULOMB_P3M_GPU:
  case COULOMB_P3M: {
#ifdef NPT
    if (q1q2) {
      double eng = p3m_add_pair_force(q1q2, d, dist2, dist, force.data());
      if (integ_switch == INTEG_METHOD_NPT_ISO)
        nptiso.p_vir[0] += eng;
    }
#else
    if (q1q2)
      p3m_add_pair_force(q1q2, d, dist2, dist, force.data());
#endif
    break;
  }
#endif
  case COULOMB_MMM1D:
    if (q1q2)
      add_mmm1d_coulomb_pair_force(q1q2, d, dist2, dist, force.data());
    break;
  case COULOMB_MMM2D:
    if (q1q2)
      add_mmm2d_coulomb_pair_force(q1q2, d, dist2, dist, force.data());
    break;
#ifdef SCAFACOS
  case COULOMB_SCAFACOS:
    if (q1q2) {
      Scafacos::add_pair_force(p1, p2, d, dist, force.data());
    }
    break;
#endif
  default:
    break;
  }

#endif /*ifdef ELECTROSTATICS */

  /***********************************************/
  /* long range magnetostatics                   */
  /***********************************************/

#ifdef DIPOLES
  /* real space magnetic dipole-dipole */
  switch (coulomb.Dmethod) {
#ifdef DP3M
  case DIPOLAR_MDLC_P3M:
  // fall trough
  case DIPOLAR_P3M: {
#ifdef NPT
    double eng = dp3m_add_pair_force(p1, p2, d, dist2, dist, force.data());
    if (integ_switch == INTEG_METHOD_NPT_ISO)
      nptiso.p_vir[0] += eng;
#else
    dp3m_add_pair_force(p1, p2, d, dist2, dist, force.data());
#endif
    break;
  }
#endif /*ifdef DP3M */
  default:
    break;
  }
#endif /* ifdef DIPOLES */

  /***********************************************/
  /* add total nonbonded forces to particle      */
  /***********************************************/

  p1->f.f += force;
  p2->f.f -= force;

  for (j = 0; j < 3; j++) {
#ifdef ROTATION
    p1->f.torque[j] += torque1[j];
    p2->f.torque[j] += torque2[j];
#endif
  }
}

/** Calculate bonded forces for one particle.
    @param p1 particle for which to calculate forces
*/

inline void add_bonded_force(Particle *p1) {
  double dx[3] = {0., 0., 0.};
  double force[3] = {0., 0., 0.};
  double force2[3] = {0., 0., 0.};
  double force3[3] = {0., 0., 0.};
#if defined(OIF_LOCAL_FORCES)
  double force4[3] = {0., 0., 0.};
#endif
#ifdef ROTATION
  double torque1[3] = {0., 0., 0.};
  double torque2[3] = {0., 0., 0.};
#endif
  Particle *p3 = nullptr, *p4 = nullptr;
  Bonded_ia_parameters *iaparams;
  int i, j, bond_broken;

  i = 0;
  while (i < p1->bl.n) {
    int type_num = p1->bl.e[i++];
    iaparams = &bonded_ia_params[type_num];
    int type = iaparams->type;
    int n_partners = iaparams->num;

    if (n_partners) {
      Particle *p2 = local_particles[p1->bl.e[i++]];
      if (!p2) {
        runtimeErrorMsg() << "bond broken between particles " << p1->p.identity;
        return;
      }

      /* fetch particle 3 eventually */
      if (n_partners >= 2) {
        p3 = local_particles[p1->bl.e[i++]];
        if (!p3) {
          runtimeErrorMsg()
              << "bond broken between particles " << p1->p.identity << ", "
              << p1->bl.e[i - 2] << " and " << p1->bl.e[i - 1]
              << " (particles are not stored on the same node)";
          return;
        }
      }

      /* fetch particle 4 eventually */
      if (n_partners >= 3) {
        p4 = local_particles[p1->bl.e[i++]];
        if (!p4) {
          runtimeErrorMsg()
              << "bond broken between particles " << p1->p.identity << ", "
              << p1->bl.e[i - 3] << ", " << p1->bl.e[i - 2] << " and "
              << p1->bl.e[i - 1] << " (particles not stored on the same node)";
          return;
        }
      }

      if (n_partners == 1) {
        /* because of the NPT pressure calculation for pair forces, we need the
           1->2 distance vector here. For many body interactions this vector is
           not needed,
           and the pressure calculation not yet clear. */
        get_mi_vector(dx, p1->r.p, p2->r.p);
      }

      switch (type) {
      case BONDED_IA_FENE:
        bond_broken = calc_fene_pair_force(p1, p2, iaparams, dx, force);
        break;
#ifdef ROTATION
      case BONDED_IA_HARMONIC_DUMBBELL:
        bond_broken =
            calc_harmonic_dumbbell_pair_force(p1, p2, iaparams, dx, force);
        break;
#endif
      case BONDED_IA_HARMONIC:
        bond_broken = calc_harmonic_pair_force(p1, p2, iaparams, dx, force);
        break;
      case BONDED_IA_QUARTIC:
        bond_broken = calc_quartic_pair_force(p1, p2, iaparams, dx, force);
        break;
      case BONDED_IA_THERMALIZED_DIST:
        bond_broken =
            calc_thermalized_bond_forces(p1, p2, iaparams, dx, force, force2);
        break;
#ifdef ELECTROSTATICS
      case BONDED_IA_BONDED_COULOMB:
        bond_broken =
            calc_bonded_coulomb_pair_force(p1, p2, iaparams, dx, force);
        break;
#endif
#ifdef P3M
      case BONDED_IA_BONDED_COULOMB_P3M_SR:
        bond_broken =
            calc_bonded_coulomb_p3m_sr_pair_force(p1, p2, iaparams, dx, force);
        break;
#endif
#ifdef MEMBRANE_COLLISION
      case BONDED_IA_OIF_OUT_DIRECTION:
        bond_broken = calc_out_direction(p1, p2, p3, p4, iaparams);
        break;
#endif
#ifdef OIF_GLOBAL_FORCES
      case BONDED_IA_OIF_GLOBAL_FORCES:
        bond_broken = 0;
        break;
#endif
#ifdef OIF_LOCAL_FORCES
      case BONDED_IA_OIF_LOCAL_FORCES:
        bond_broken = calc_oif_local(p1, p2, p3, p4, iaparams, force, force2,
                                     force3, force4);
        break;
#endif
// IMMERSED_BOUNDARY
#ifdef IMMERSED_BOUNDARY
      /*      case BONDED_IA_IBM_WALL_REPULSION:
              IBM_WallRepulsion_CalcForce(p1, iaparams);
              bond_broken = 0;
              // These may be added later on, but we set them to zero because
         the force has already been added in IBM_WallRepulsion_CalcForce
              force[0] = force2[0] = force3[0] = 0;
              force[1] = force2[1] = force3[1] = 0;
              force[2] = force2[2] = force3[2] = 0;
              break;*/
      case BONDED_IA_IBM_TRIEL:
        bond_broken = IBM_Triel_CalcForce(p1, p2, p3, iaparams);
        // These may be added later on, but we set them to zero because the
        // force has already been added in IBM_Triel_CalcForce
        force[0] = force2[0] = force3[0] = 0;
        force[1] = force2[1] = force3[1] = 0;
        force[2] = force2[2] = force3[2] = 0;
        break;
      case BONDED_IA_IBM_VOLUME_CONSERVATION:
        bond_broken = 0;
        // Don't do anything here. We calculate and add the global volume forces
        // in IBM_VolumeConservation. They cannot be calculated on a per-bond
        // basis
        force[0] = force2[0] = force3[0] = 0;
        force[1] = force2[1] = force3[1] = 0;
        force[2] = force2[2] = force3[2] = 0;
        break;
      case BONDED_IA_IBM_TRIBEND: {
        // First build neighbor list. This includes all nodes around the central
        // node.
        const int numNeighbors = iaparams->num;
        Particle **neighbors = new Particle *[numNeighbors];
        // Three are already there
        neighbors[0] = p2;
        neighbors[1] = p3;
        neighbors[2] = p4;
        // Get rest
        for (int j = 3; j < numNeighbors; j++)
          neighbors[j] = local_particles[p1->bl.e[i++]];

        IBM_Tribend_CalcForce(p1, numNeighbors, neighbors, *iaparams);
        bond_broken = 0;

        // Clean up
        delete[] neighbors;

        // These may be added later on, but we set them to zero because the
        // force has
        force[0] = force2[0] = force3[0] = 0;
        force[1] = force2[1] = force3[1] = 0;
        force[2] = force2[2] = force3[2] = 0;
        break;
      }
#endif

#ifdef LENNARD_JONES
      case BONDED_IA_SUBT_LJ:
        bond_broken = calc_subt_lj_pair_force(p1, p2, iaparams, dx, force);
        break;
#endif
#ifdef BOND_ANGLE
      case BONDED_IA_ANGLE_HARMONIC:
        bond_broken =
            calc_angle_harmonic_force(p1, p2, p3, iaparams, force, force2);
        break;
      case BONDED_IA_ANGLE_COSINE:
        bond_broken =
            calc_angle_cosine_force(p1, p2, p3, iaparams, force, force2);
        break;
      case BONDED_IA_ANGLE_COSSQUARE:
        bond_broken =
            calc_angle_cossquare_force(p1, p2, p3, iaparams, force, force2);
        break;
#endif
#ifdef BOND_ANGLEDIST
      case BONDED_IA_ANGLEDIST:
        bond_broken = calc_angledist_force(p1, p2, p3, iaparams, force, force2);
        break;
#endif
      case BONDED_IA_DIHEDRAL:
        bond_broken = calc_dihedral_force(p1, p2, p3, p4, iaparams, force,
                                          force2, force3);
        break;
#ifdef BOND_CONSTRAINT
      case BONDED_IA_RIGID_BOND:
        // add_rigid_bond_pair_force(p1,p2, iaparams, force, force2);
        bond_broken = 0;
        force[0] = force[1] = force[2] = 0.0;
        break;
#endif
#ifdef TABULATED
      case BONDED_IA_TABULATED:
        switch (iaparams->p.tab.type) {
        case TAB_BOND_LENGTH:
          bond_broken = calc_tab_bond_force(p1, p2, iaparams, dx, force);
          break;
        case TAB_BOND_ANGLE:
          bond_broken =
              calc_tab_angle_force(p1, p2, p3, iaparams, force, force2);
          break;
        case TAB_BOND_DIHEDRAL:
          bond_broken = calc_tab_dihedral_force(p1, p2, p3, p4, iaparams, force,
                                                force2, force3);
          break;
        default:
          runtimeErrorMsg() << "add_bonded_force: tabulated bond type of atom "
                            << p1->p.identity << " unknown\n";
          return;
        }
        break;
#endif
#ifdef UMBRELLA
      case BONDED_IA_UMBRELLA:
        bond_broken = calc_umbrella_pair_force(p1, p2, iaparams, dx, force);
        break;
#endif
      case BONDED_IA_VIRTUAL_BOND:
        bond_broken = 0;
        force[0] = force[1] = force[2] = 0.0;
        break;
      default:
        runtimeErrorMsg() << "add_bonded_force: tabulated bond type of atom "
                          << p1->p.identity << " unknown\n";
        return;
      }

      switch (n_partners) {
      case 1:
        if (bond_broken) {
          runtimeErrorMsg()
              << "bond broken between particles " << p1->p.identity << " and "
              << p2->p.identity << ". Distance vector: " << dx[0] << " "
              << dx[1] << " " << dx[2];
          continue;
        }

        for (j = 0; j < 3; j++) {
          switch (type) {
          case BONDED_IA_THERMALIZED_DIST:
            p1->f.f[j] += force[j];
            p2->f.f[j] += force2[j];
            break;
          default:
            p1->f.f[j] += force[j];
            p2->f.f[j] -= force[j];
#ifdef ROTATION
            p1->f.torque[j] += torque1[j];
            p2->f.torque[j] += torque2[j];
#endif
          }

#ifdef NPT
          if (integ_switch == INTEG_METHOD_NPT_ISO)
            nptiso.p_vir[j] += force[j] * dx[j];
#endif
        }
        break;
      case 2:
        if (bond_broken) {
          runtimeErrorMsg()
              << "bond broken between particles " << p1->p.identity << ", "
              << p2->p.identity << " and " << p3->p.identity;
          continue;
        }

        for (j = 0; j < 3; j++) {
          switch (type) {
#ifdef OIF_GLOBAL_FORCES
          case BONDED_IA_OIF_GLOBAL_FORCES:
            break;
#endif
          default:
            p1->f.f[j] += force[j];
            p2->f.f[j] += force2[j];
            p3->f.f[j] -= (force[j] + force2[j]);
          }
        }
        break;
      case 3:
        if (bond_broken) {
          runtimeErrorMsg() << "bond broken between particles "
                            << p1->p.identity << ", " << p2->p.identity << ", "
                            << p3->p.identity << " and " << p4->p.identity;
          continue;
        }

        switch (type) {
        case BONDED_IA_DIHEDRAL:
          for (j = 0; j < 3; j++) {
            p1->f.f[j] += force[j];
            p2->f.f[j] += force2[j];
            p3->f.f[j] += force3[j];
            p4->f.f[j] -= force[j] + force2[j] + force3[j];
          }
          break;

#ifdef OIF_LOCAL_FORCES
        case BONDED_IA_OIF_LOCAL_FORCES:
          for (j = 0; j < 3; j++) {
            p1->f.f[j] += force2[j];
            p2->f.f[j] += force[j];
            p3->f.f[j] += force3[j];
            p4->f.f[j] += force4[j];
          }
          break;
#endif
        }
        break;
      }
    }
  }
}

inline void check_particle_force(Particle *part) {
  for (int i = 0; i < 3; i++) {
    if (std::isnan(part->f.f[i])) {
      runtimeErrorMsg() << "force on particle " << part->p.identity
                        << " was NAN.";
    }
  }

#ifdef ROTATION
  for (int i = 0; i < 3; i++) {
    if (std::isnan(part->f.torque[i])) {
      runtimeErrorMsg() << "torque on particle " << part->p.identity
                        << " was NAN.";
    }
  }
#endif
}

inline void add_single_particle_force(Particle *p) {
  add_bonded_force(p);
}

#endif
