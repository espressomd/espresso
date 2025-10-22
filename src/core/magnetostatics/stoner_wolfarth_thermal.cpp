/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#include "config/config.hpp"

#ifdef DIPOLES
#define TWO_M_PI 2 * M_PI

#include "magnetostatics/dipolar_direct_sum.hpp"

#include "cells.hpp"
#include "communication.hpp"
#include "constraints.hpp"
#include "constraints/HomogeneousMagneticField.hpp"
#include "errorhandling.hpp"
#include "grid.hpp"

#include <utils/cartesian_product.hpp>
#include <utils/constants.hpp>
#include <utils/math/sqr.hpp>
#include <utils/math/vec_rotate.hpp>
#include <utils/mpi/iall_gatherv.hpp>

#include <boost/mpi/collectives.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/range/counting_range.hpp>

#include <nlopt.hpp>

#include "event.hpp"
#include "magnetostatics/stoner_wolfarth_thermal.hpp"
#include "rotation.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iterator>
#include <mpi.h>
#include <random>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace {
double phi_objective(unsigned n, const double *x, double *grad,
                     void *my_func_data) {
  double phi = x[0];
  double *params = (double *)my_func_data;
  double theta = params[0];
  double h = params[1];
  if (grad) {
    grad[0] = std::sin(2 * (phi - theta)) + 2 * h * std::sin(phi);
  }
  return -0.5 - 0.5 * std::cos(2 * (phi - theta)) - 2 * h * std::cos(phi);
}

/*
SW energy minimisation step with MC dip_mom flip step. SW energy normalised by
the anisotropy field H_k (Hkinv stored on part to avoid division, h is the
reduced field  due to the normalisation)
*/
double funct(double theta, double h, double phi0, double kT_KVm_inv,
             double tau0_inv, double dt, std::mt19937 &rng_generator) {

  std::uniform_real_distribution<double> distribution(0.0, 1.0);
  double eps_phi = 1e-3;
  double h_crit = std::pow(std::pow(std::sin(theta), 2.0 / 3) +
                               std::pow(std::cos(theta), 2.0 / 3),
                           -3.0 / 2);
  nlopt::opt opt(nlopt::LD_MMA, 1);
  double params[] = {theta, h};

  opt.set_min_objective(phi_objective, &params);
  opt.set_ftol_rel(
      1e-15); // Set the relative tolerance for the objective function value
  opt.set_ftol_abs(
      1e-15); // Set the relative tolerance for the objective function value
  std::vector<double> x(1);

  x[0] = phi0 + eps_phi; /* make initial guess from previos position plus
                                   an arbitrary perturbation*/
  double min1; /* this is the actuall value of the energy from minimiser */
  opt.optimize(x, min1);
  double phi_min1 = fmod(x[0], TWO_M_PI);
  double sol = phi_min1;
  if (fabs(h) < h_crit) {
    opt.set_max_objective(phi_objective, &params);
    x[0] = phi0 + eps_phi;
    double max1;
    opt.optimize(x, max1);
    double phi_max1 = fmod(x[0], TWO_M_PI);
    x[0] = fmod(phi_max1 + M_PI, 2 * M_PI);
    double max2;
    opt.optimize(x, max2);

    double b1 = std::abs(max1 - min1) * kT_KVm_inv;
    double b2 = std::abs(max2 - min1) * kT_KVm_inv;
    double b_min = (b1 < b2) ? b1 : b2;
    double tau_inv = tau0_inv * exp(-b_min);
    double p12 = 1. - exp(-dt * tau_inv);

    if (distribution(rng_generator) < p12) {
      opt.set_min_objective(phi_objective, &params);
      x[0] = fmod(phi_min1 + M_PI + eps_phi, 2 * M_PI);
      /*try to find another minimimum from the other side*/
      double min2;
      opt.optimize(x, min2);

      double phi_min2 = fmod(x[0], TWO_M_PI);
      sol = phi_min2;
    }
  }
  return fmod(sol + TWO_M_PI, TWO_M_PI);
}

double funct_tans(double sol, double theta, double h, double kT_KVm_inv,
                  double tau_trans_inv, double dt,
                  std::mt19937 &rng_generator) {

  std::uniform_real_distribution<double> ang_distribution(-M_PI_2, M_PI_2);
  std::uniform_real_distribution<double> distribution(0, 1);

  double shift = ang_distribution(rng_generator);
  double attempot_ang = sol + shift;
  double params[] = {theta, h};
  double b_min = std::abs(phi_objective(1, &sol, nullptr, &params) -
                          phi_objective(1, &attempot_ang, nullptr, &params)) *
                 kT_KVm_inv;
  double tau_inv = tau_trans_inv * b_min;
  double p12 = 1. - exp(-dt * tau_inv);
  // std::cout << p12 << "\n";s
  if (distribution(rng_generator) > p12) {
    shift = 0;
  }
  return shift;
}

} // namespace

void stoner_wolfarth_main(ParticleRange const &particles,
                          std::mt19937 &rng_generator) {
  /* collect particle data */
  std::vector<Particle *> local_real_particles;
  std::vector<Particle *> local_virt_particles;
  local_real_particles.reserve(particles.size());
  local_virt_particles.reserve(particles.size());
  for (auto &p : particles) {
    if (p.sw_real() == 1) {
      local_real_particles.emplace_back(&p);
    } else if (p.sw_virt() == 1) {
      local_virt_particles.emplace_back(&p);
    }
  }
  // must assert that there is an equal number of sw_reals and sw_virts
  Utils::Vector3d cntrl = {0., 0., 0.};
  Utils::Vector3d ext_fld = {0., 0., 0.};
  /* collect HomogeneousMagneticFields if active */
  for (auto const &constraint : ::Constraints::constraints) {
    auto ptr = dynamic_cast<::Constraints::HomogeneousMagneticField *const>(
        &*constraint);
    if (ptr != nullptr) {
      ext_fld += ptr->H();
    }
  }
  if (ext_fld != cntrl) {
    auto p = local_virt_particles.begin();
    for (auto pi = local_real_particles.begin();
         pi != local_real_particles.end(); ++pi, ++p) {
      Utils::Vector3d ext_fld_dpl = {0., 0., 0.};
      ext_fld_dpl = ext_fld + (*p)->dip_fld();
      double h = ext_fld_dpl.norm() * (*p)->Hkinv();
      auto e_h = ext_fld_dpl.normalized();
      // calc_director() result already normalised
      Utils::Vector3d e_k = (*pi)->calc_director();
      double theta = std::acos(e_h * e_k);
      if (theta > M_PI_2) {
        theta = M_PI - theta;
        h = -h;
        e_h = -e_h;
      }
      auto rot_axis =
          vector_product(vector_product(e_h, e_k), e_h).normalized();
      auto phi = funct(theta, h, (*pi)->phi0(), (*pi)->kT_KVm_inv(),
                       (*pi)->tau0_inv(), (*pi)->dt_incr(), rng_generator);
      auto shift =
          funct_tans(phi, theta, h, (*pi)->kT_KVm_inv(), (*pi)->tau_trans_inv(),
                     (*pi)->dt_incr(), rng_generator);
      auto mom = e_h * std::cos(phi) + rot_axis * std::sin(phi);
      mom = Utils::vec_rotate(vector_product(e_h, rot_axis), shift, mom);
      (*pi)->phi0() = phi + shift;
      auto const [quat, dipm] = convert_dip_to_quat(mom * (*p)->sat_mag());
      (*p)->dipm() = dipm;
      (*p)->quat() = quat;
    }
    // on_dipoles_change();
  } else {
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    auto p = local_virt_particles.begin();
    for (auto pi = local_real_particles.begin();
         pi != local_real_particles.end(); ++pi, ++p) {
      Utils::Vector3d e_k = (*pi)->calc_director();
      double tau_inv = (*pi)->tau0_inv() * exp(-(*pi)->kT_KVm_inv());
      double p12 = 1. - exp(-(*pi)->dt_incr() * tau_inv);
      if (distribution(rng_generator) < p12) {
        if ((*pi)->phi0() == 0) {
          auto const [quat, dipm] = convert_dip_to_quat((*p)->sat_mag() * -e_k);
          (*pi)->phi0() = M_PI;
          (*p)->dipm() = dipm;
          (*p)->quat() = quat;
        } else if ((*pi)->phi0() == M_PI) {
          auto const [quat, dipm] = convert_dip_to_quat((*p)->sat_mag() * e_k);
          (*pi)->phi0() = 0;
          (*p)->dipm() = dipm;
          (*p)->quat() = quat;
        } else {
          double diff_0 = std::abs((*pi)->phi0() - 0);
          double diff_PI = std::abs((*pi)->phi0() - M_PI);
          // Compare the differences and determine the closer angle
          if (diff_0 < diff_PI) {
            auto const [quat, dipm] =
                convert_dip_to_quat((*p)->sat_mag() * -e_k);
            (*pi)->phi0() = M_PI;
            (*p)->dipm() = dipm;
            (*p)->quat() = quat;
          } else {
            auto const [quat, dipm] =
                convert_dip_to_quat((*p)->sat_mag() * e_k);
            (*pi)->phi0() = 0;
            (*p)->dipm() = dipm;
            (*p)->quat() = quat;
          }
        }
      } else {
        if ((*pi)->phi0() == 0) {
          auto const [quat, dipm] = convert_dip_to_quat((*p)->sat_mag() * e_k);
          (*p)->dipm() = dipm;
          (*p)->quat() = quat;
        } else if ((*pi)->phi0() == M_PI) {
          auto const [quat, dipm] = convert_dip_to_quat((*p)->sat_mag() * -e_k);
          (*p)->dipm() = dipm;
          (*p)->quat() = quat;
        } else {
          double diff_0 = std::abs((*pi)->phi0() - 0);
          double diff_PI = std::abs((*pi)->phi0() - M_PI);
          // Compare the differences and determine the closer angle
          if (diff_0 < diff_PI) {
            auto const [quat, dipm] =
                convert_dip_to_quat((*p)->sat_mag() * e_k);
            (*pi)->phi0() = 0;
            (*p)->dipm() = dipm;
            (*p)->quat() = quat;
          } else {
            auto const [quat, dipm] =
                convert_dip_to_quat((*p)->sat_mag() * -e_k);
            (*pi)->phi0() = M_PI;
            (*p)->dipm() = dipm;
            (*p)->quat() = quat;
          }
        }
      }
    }
  }
  // this call might be necessart when using p3m! significant overhead
  // on_dipoles_change();
}
#endif