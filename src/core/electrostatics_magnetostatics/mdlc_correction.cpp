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
 */

#include "electrostatics_magnetostatics/mdlc_correction.hpp"

#ifdef DP3M
#include "Particle.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "electrostatics_magnetostatics/dipole.hpp"
#include "electrostatics_magnetostatics/p3m-dipolar.hpp"
#include "errorhandling.hpp"
#include "global.hpp"
#include "grid.hpp"
#include "particle_data.hpp"

DLC_struct dlc_params = {1e100, 0, 0, 0, 0};

static double mu_max;

void calc_mu_max() {
  auto local_particles = cell_structure.local_particles();
  mu_max = std::accumulate(
      local_particles.begin(), local_particles.end(), 0.0,
      [](double mu, Particle const &p) { return std::max(mu, p.p.dipm); });

  MPI_Allreduce(MPI_IN_PLACE, &mu_max, 1, MPI_DOUBLE, MPI_MAX, comm_cart);
}

inline double g1_DLC_dip(double g, double x) {
  auto const c = g / x;
  auto const cc2 = c * c;
  auto const x3 = x * x * x;
  auto const a = g * g * g / x + 1.5 * cc2 + 1.5 * g / x3 + 0.75 / (x3 * x);
  return a;
}

inline double g2_DLC_dip(double g, double x) {
  auto const x2 = x * x;
  auto const a = g * g / x + 2.0 * g / x2 + 2.0 / (x2 * x);
  return a;
}

/** Compute Mx, My, Mz and Mtotal */
double slab_dip_count_mu(double *mt, double *mx, double *my,
                         const ParticleRange &particles) {
  Utils::Vector3d node_sums{};
  Utils::Vector3d tot_sums{};

  for (auto const &p : particles) {
    if (p.p.dipm != 0.0) {
      node_sums += p.calc_dip();
    }
  }

  MPI_Allreduce(node_sums.data(), tot_sums.data(), 3, MPI_DOUBLE, MPI_SUM,
                comm_cart);

  auto const M = tot_sums.norm();
  auto const Mz = tot_sums[2];
  auto const Mx = tot_sums[0];
  auto const My = tot_sums[1];

  *mt = M;
  *mx = Mx;
  *my = My;

  return Mz;
}

/** Compute the dipolar DLC corrections for forces and torques.
 *  %Algorithm implemented accordingly to @cite brodka04a.
 */
double get_DLC_dipolar(int kcut, std::vector<Utils::Vector3d> &fs,
                       std::vector<Utils::Vector3d> &ts,
                       const ParticleRange &particles) {
  auto const n_local_particles = particles.size();

  std::vector<double> ReSjp(n_local_particles), ReSjm(n_local_particles);
  std::vector<double> ImSjp(n_local_particles), ImSjm(n_local_particles);
  std::vector<double> ReGrad_Mup(n_local_particles),
      ImGrad_Mup(n_local_particles);
  std::vector<double> ReGrad_Mum(n_local_particles),
      ImGrad_Mum(n_local_particles);
  double s1, s2, s3, s4;
  double s1z, s2z, s3z, s4z;
  double ss;

  auto const facux = 2.0 * M_PI / box_geo.length()[0];
  auto const facuy = 2.0 * M_PI / box_geo.length()[1];
  double energy = 0.0;

  for (int ix = -kcut; ix <= +kcut; ix++) {
    for (int iy = -kcut; iy <= +kcut; iy++) {
      if (!(ix == 0 && iy == 0)) {
        auto const gx = static_cast<double>(ix) * facux;
        auto const gy = static_cast<double>(iy) * facuy;

        auto const gr = sqrt(gx * gx + gy * gy);

        // We assume short slab direction is z direction
        auto const fa1 = 1. / (gr * (exp(gr * box_geo.length()[2]) - 1.0));

        // ... Compute S+,(S+)*,S-,(S-)*, and Spj,Smj for the current g

        double S[4] = {0.0, 0.0, 0.0, 0.0}; // S of Brodka method, or is S[4] =
                                            // {Re(S+), Im(S+), Re(S-), Im(S-)}
        int ip = 0;

        for (auto const &p : particles) {
          if (p.p.dipm > 0) {
            Utils::Vector3d const dip = p.calc_dip();

            auto const a = gx * dip[0] + gy * dip[1];
            auto const b = gr * dip[2];
            auto const er = gx * p.r.p[0] + gy * p.r.p[1];
            auto const ez = gr * p.r.p[2];
            auto const c = cos(er);
            auto const d = sin(er);
            auto const f = exp(ez);

            ReSjp[ip] = (b * c - a * d) * f;
            ImSjp[ip] = (c * a + b * d) * f;
            ReSjm[ip] = (-b * c - a * d) / f;
            ImSjm[ip] = (c * a - b * d) / f;
            ReGrad_Mup[ip] = c * f;
            ReGrad_Mum[ip] = c / f;
            ImGrad_Mup[ip] = d * f;
            ImGrad_Mum[ip] = d / f;

            S[0] += ReSjp[ip];
            S[1] += ImSjp[ip];
            S[2] += ReSjm[ip];
            S[3] += ImSjm[ip];
          }
          ip++;
        }

        MPI_Allreduce(MPI_IN_PLACE, S, 4, MPI_DOUBLE, MPI_SUM, comm_cart);

        // We compute the contribution to the energy ............

        // s2=(ReSm*ReSp+ImSm*ImSp); s2=s1!!!

        energy += fa1 * ((S[0] * S[2] + S[1] * S[3]) * 2.0);

        // ... Now we can compute the contributions to E,Fj,Ej for the current
        // g-value
        ip = 0;
        for (auto &p : particles) {
          if (p.p.dipm > 0) {
            // We compute the contributions to the forces ............

            s1 = -(-ReSjp[ip] * S[3] + ImSjp[ip] * S[2]);
            s2 = +(ReSjm[ip] * S[1] - ImSjm[ip] * S[0]);
            s3 = -(-ReSjm[ip] * S[1] + ImSjm[ip] * S[0]);
            s4 = +(ReSjp[ip] * S[3] - ImSjp[ip] * S[2]);

            s1z = +(ReSjp[ip] * S[2] + ImSjp[ip] * S[3]);
            s2z = -(ReSjm[ip] * S[0] + ImSjm[ip] * S[1]);
            s3z = -(ReSjm[ip] * S[0] + ImSjm[ip] * S[1]);
            s4z = +(ReSjp[ip] * S[2] + ImSjp[ip] * S[3]);

            ss = s1 + s2 + s3 + s4;
            fs[ip][0] += fa1 * gx * ss;
            fs[ip][1] += fa1 * gy * ss;
            fs[ip][2] += fa1 * gr * (s1z + s2z + s3z + s4z);

            // We compute the contributions to the electrical field
            // ............

            s1 = -(-ReGrad_Mup[ip] * S[3] + ImGrad_Mup[ip] * S[2]);
            s2 = +(ReGrad_Mum[ip] * S[1] - ImGrad_Mum[ip] * S[0]);
            s3 = -(-ReGrad_Mum[ip] * S[1] + ImGrad_Mum[ip] * S[0]);
            s4 = +(ReGrad_Mup[ip] * S[3] - ImGrad_Mup[ip] * S[2]);

            s1z = +(ReGrad_Mup[ip] * S[2] + ImGrad_Mup[ip] * S[3]);
            s2z = -(ReGrad_Mum[ip] * S[0] + ImGrad_Mum[ip] * S[1]);
            s3z = -(ReGrad_Mum[ip] * S[0] + ImGrad_Mum[ip] * S[1]);
            s4z = +(ReGrad_Mup[ip] * S[2] + ImGrad_Mup[ip] * S[3]);

            ss = s1 + s2 + s3 + s4;
            ts[ip][0] += fa1 * gx * ss;
            ts[ip][1] += fa1 * gy * ss;
            ts[ip][2] += fa1 * gr * (s1z + s2z + s3z + s4z);
          } // if dipm>0 ....
          ip++;
        } // loop j
      }   // end of if(ii> ...
    }
  } // end of loops for gx,gy

  // Convert from the corrections to the Electrical field to the corrections
  // for the torques ....

  int ip = 0;
  for (auto const &p : particles) {
    if (p.p.dipm > 0) {
      ts[ip] = vector_product(p.calc_dip(), ts[ip]);
    }
    ip++;
  }

  // Multiply by the factors we have left during the loops

  auto const piarea = M_PI / (box_geo.length()[0] * box_geo.length()[1]);

  for (int j = 0; j < n_local_particles; j++) {
    fs[j] *= piarea;
    ts[j] *= piarea;
  }

  energy *= (-piarea);

  return energy;
}

/** Compute the dipolar DLC corrections
 *  %Algorithm implemented accordingly to @cite brodka04a.
 */
double get_DLC_energy_dipolar(int kcut, const ParticleRange &particles) {
  auto const facux = 2.0 * M_PI / box_geo.length()[0];
  auto const facuy = 2.0 * M_PI / box_geo.length()[1];

  double energy = 0.0;
  for (int ix = -kcut; ix <= +kcut; ix++) {
    for (int iy = -kcut; iy <= +kcut; iy++) {

      if (!(ix == 0 && iy == 0)) {
        auto const gx = static_cast<double>(ix) * facux;
        auto const gy = static_cast<double>(iy) * facuy;
        auto const gr = sqrt(gx * gx + gy * gy);
        // We assume short slab direction is z direction
        auto const fa1 = 1. / (gr * (exp(gr * box_geo.length()[2]) - 1.0));

        // ... Compute S+,(S+)*,S-,(S-)*, and Spj,Smj for the current g

        double S[4] = {0.0, 0.0, 0.0, 0.0};
        int ip = 0;
        for (auto const &p : particles) {
          if (p.p.dipm > 0) {
            const Utils::Vector3d dip = p.calc_dip();

            auto const a = gx * dip[0] + gy * dip[1];
            {
              auto const b = gr * dip[2];
              auto const er = gx * p.r.p[0] + gy * p.r.p[1];
              auto const ez = gr * p.r.p[2];
              auto const c = cos(er);
              auto const d = sin(er);
              auto const f = exp(ez);

              S[0] += (b * c - a * d) * f;
              S[1] += (c * a + b * d) * f;
              S[2] += (-b * c - a * d) / f;
              S[3] += (c * a - b * d) / f;
            }
          }
          ip++;
        }

        double global_S[4];
        MPI_Reduce(S, global_S, 4, MPI_DOUBLE, MPI_SUM, 0, comm_cart);

        // We compute the contribution to the energy ............
        auto const s1 = global_S[0] * global_S[2] + global_S[1] * global_S[3];
        // s2=(ReSm*ReSp+ImSm*ImSp); s2=s1!!!

        energy += fa1 * (s1 * 2.0);

      } // end of if(...
    }
  } // end of loops for gx,gy

  // Multiply by the factors we have left during the loops

  auto const piarea = M_PI / (box_geo.length()[0] * box_geo.length()[1]);
  energy *= (-piarea);
  return (this_node == 0) ? energy : 0.0;
}

/** Compute and add the terms needed to correct the 3D dipolar
 *  methods when we have a slab geometry
 */
void add_mdlc_force_corrections(const ParticleRange &particles) {
  auto const volume = box_geo.volume();

  // --- Create arrays that should contain the corrections to
  //     the forces and torques, and set them to zero.
  std::vector<Utils::Vector3d> dip_DLC_f(particles.size());
  std::vector<Utils::Vector3d> dip_DLC_t(particles.size());

  //---- Compute the corrections ----------------------------------

  // First the DLC correction
  get_DLC_dipolar(static_cast<int>(std::round(dlc_params.far_cut)), dip_DLC_f,
                  dip_DLC_t, particles);

  // Now we compute the correction like Yeh and Klapp to take into account
  // the fact that you are using a 3D PBC method which uses spherical
  // summation instead of slab-wise summation.
  // Slab-wise summation is the one required to apply DLC correction.
  // This correction is often called SDC = Shape Dependent Correction.
  // See @cite brodka04a.

  double mx = 0.0, my = 0.0, mtot = 0.0;
  auto const mz = slab_dip_count_mu(&mtot, &mx, &my, particles);

  // --- Transfer the computed corrections to the Forces, Energy and torques
  //     of the particles

  int ip = 0;
  for (auto &p : particles) {
    if ((p.p.dipm) != 0.0) {
      // SDC correction term is zero for the forces
      p.f.f += dipole.prefactor * dip_DLC_f[ip];

#if defined(ROTATION) && defined(DP3M)
      auto const dip = p.calc_dip();
      auto const correc = 4. * M_PI / volume;
      Utils::Vector3d d;
      // in the Next lines: the second term (correc*...) is the SDC
      // correction for the torques
      if (dp3m.params.epsilon == P3M_EPSILON_METALLIC) {
        d = {0.0, 0.0, -correc * mz};
      } else {
        auto const correps = correc / (2.0 * dp3m.params.epsilon + 1.0);
        d = {correps * mx, correps * my,
             correc * mz * (-1.0 + 1. / (2.0 * dp3m.params.epsilon + 1.0))};
      }
      p.f.torque += dipole.prefactor * (dip_DLC_t[ip] + vector_product(dip, d));
#endif
    }
    ip++;
  }
}

/** Compute and add the terms needed to correct the energy of
 *  3D dipolar methods when we have a slab geometry
 */
double add_mdlc_energy_corrections(const ParticleRange &particles) {

  auto const volume = box_geo.volume();

  //---- Compute the corrections ----------------------------------

  // First the DLC correction
  double dip_DLC_energy =
      dipole.prefactor *
      get_DLC_energy_dipolar(static_cast<int>(std::round(dlc_params.far_cut)),
                             particles);

  // Now we compute the correction like Yeh and Klapp to take into account
  // the fact that you are using a 3D PBC method which uses spherical
  // summation instead of slab-wise summation.
  // Slab-wise summation is the one required to apply DLC correction.
  // This correction is often called SDC = Shape Dependent Correction.
  // See @cite brodka04a.

  double mx = 0.0, my = 0.0, mtot = 0.0;
  auto const mz = slab_dip_count_mu(&mtot, &mx, &my, particles);
  auto const mz2 = mz * mz;

  if (this_node == 0) {
#ifdef DP3M
    if (dipole.method == DIPOLAR_MDLC_P3M) {
      if (dp3m.params.epsilon == P3M_EPSILON_METALLIC) {
        dip_DLC_energy += dipole.prefactor * 2. * M_PI / volume * mz2;
      } else {
        dip_DLC_energy +=
            dipole.prefactor * 2. * M_PI / volume *
            (mz2 - mtot * mtot / (2.0 * dp3m.params.epsilon + 1.0));
      }
    } else
#endif
    {
      dip_DLC_energy += dipole.prefactor * 2. * M_PI / volume * mz2;
      fprintf(stderr, "You are not using the P3M method, therefore "
                      "dp3m.params.epsilon unknown, I assume metallic borders "
                      "\n");
    }

    return dip_DLC_energy;
  }
  return 0.0;
}

/** Compute the cut-off in the DLC dipolar part to get a certain accuracy.
 *  We assume particles to have all the same value of the dipolar momentum
 *  modulus (@ref mu_max). @ref mu_max is taken as the largest value of mu
 *  inside the system. If we assume the gap has a width @c gap_size (within
 *  which there is no particles): <tt>Lz = h + gap_size</tt>
 *
 *  BE CAREFUL:
 *  1. We assume the short distance for the slab to be in the *z*-direction
 *  2. You must also tune the other 3D method to the same accuracy, otherwise
 *     it makes no sense to have an accurate result for DLC-dipolar.
 */
int mdlc_tune(double error) {
  mpi_bcast_max_mu(); /* we take the maximum dipole in the system, to be sure
                         that the errors in the other case
                         will be equal or less than for this one */

  const double n = get_n_part();
  auto const lx = box_geo.length()[0];
  auto const ly = box_geo.length()[1];
  auto const lz = box_geo.length()[2];
  auto const a = lx * ly;
  auto const h = dlc_params.h;
  if (h < 0)
    return ES_ERROR;

  if (h > lz) {
    fprintf(stderr,
            "tune DLC dipolar: Slab is larger than the box size !!! \n");
    errexit();
  }

  if (fabs(box_geo.length()[0] - box_geo.length()[1]) > 0.001) {
    fprintf(stderr, "tune DLC dipolar: box size in x direction is different "
                    "from y direction !!! \n");
    fprintf(stderr, "The tuning formula requires both to be equal. \n");
    errexit();
  }

  bool flag = false;
  int kc;
  int const limitkc = 200;
  for (kc = 1; kc < limitkc; kc++) {
    auto const gc = kc * 2.0 * Utils::pi() / lx;
    auto const fa0 = sqrt(9.0 * exp(+2. * gc * h) * g1_DLC_dip(gc, lz - h) +
                          22.0 * g1_DLC_dip(gc, lz) +
                          9.0 * exp(-2.0 * gc * h) * g1_DLC_dip(gc, lz + h));
    auto const fa1 = 0.5 * sqrt(Utils::pi() / (2.0 * a)) * fa0;
    auto const fa2 = g2_DLC_dip(gc, lz);
    auto const de =
        n * (mu_max * mu_max) / (4.0 * (exp(gc * lz) - 1.0)) * (fa1 + fa2);
    if (de < error) {
      flag = true;
      break;
    }
  }

  if (!flag) {
    fprintf(stderr, "tune DLC dipolar: Sorry, unable to find a proper cut-off "
                    "for such system and accuracy.\n");
    fprintf(stderr, "Try modifying the variable limitkc in the c-code: "
                    "mdlc_correction.cpp  ... \n");
    return ES_ERROR;
  }

  dlc_params.far_cut = kc;

  return ES_OK;
}

int mdlc_sanity_checks() {
  if (!box_geo.periodic(0) || !box_geo.periodic(1) || !box_geo.periodic(2)) {
    runtimeErrorMsg() << "mdlc requires periodicity 1 1 1";
    return 1;
  }

  // It will be desirable to have a  checking function that check that the
  // slab
  // geometry is such that
  // the short direction is along the z component.

  return 0;
}

int mdlc_set_params(double maxPWerror, double gap_size, double far_cut) {

  dlc_params.maxPWerror = maxPWerror;
  dlc_params.gap_size = gap_size;
  dlc_params.h = box_geo.length()[2] - gap_size;

  if (Dipole::set_mesh()) {
    // if Dipole::set_mesh fails
    return ES_ERROR;
  }

  dlc_params.far_cut = far_cut;
  if (far_cut != -1) {
    dlc_params.far_calculated = 0;
  } else {
    dlc_params.far_calculated = 1;
    if (mdlc_tune(dlc_params.maxPWerror) == ES_ERROR) {
      runtimeErrorMsg() << "mdlc tuning failed, gap size too small";
    }
  }
  mpi_bcast_coulomb_params();

  return ES_OK;
}

#endif // DP3M
