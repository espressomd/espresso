/*
 * Copyright (C) 2014-2022 The ESPResSo project
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

/**
 * @file
 * MMM1D algorithm for long-range %Coulomb interactions on the GPU.
 * Implementation of the MMM1D method for the calculation of the electrostatic
 * interaction in one-dimensionally periodic systems. For details on the
 * method see MMM in general. The MMM1D method works only with the N-squared
 * cell system since neither the near nor far formula can be decomposed.
 */

#ifndef ESPRESSO_SRC_CORE_ELECTROSTATICS_MMM1D_GPU_HPP
#define ESPRESSO_SRC_CORE_ELECTROSTATICS_MMM1D_GPU_HPP

#include "config/config.hpp"

#ifdef MMM1D_GPU

#include "electrostatics/actor.hpp"

class CoulombMMM1DGpu : public Coulomb::Actor<CoulombMMM1DGpu> {
public:
  double maxPWerror;
  double far_switch_radius;
  double far_switch_radius_sq;
  int bessel_cutoff;

  CoulombMMM1DGpu(double prefactor, double maxPWerror, double far_switch_radius,
                  int bessel_cutoff);
  ~CoulombMMM1DGpu();

  // interface methods
  void add_long_range_forces();
  void add_long_range_energy();

  void on_activation() {
    sanity_checks();
    tune();
  }
  void on_boxl_change() { setup(); }
  void on_node_grid_change() const {}
  void on_periodicity_change() const { sanity_checks_periodicity(); }
  void on_cell_structure_change() { sanity_checks_cell_structure(); }
  void init() const {}

  void sanity_checks() const {
    sanity_checks_periodicity();
    sanity_checks_cell_structure();
    sanity_checks_charge_neutrality();
  }

  void tune();
  bool is_tuned() const { return m_is_tuned; }

private:
  bool m_is_tuned;

  // the box length currently set on the GPU
  // Needed to make sure it hasn't been modified after inter coulomb was used.
  float host_boxz = 0.f;
  // the number of particles we had during the last run. Needed to check if we
  // have to realloc dev_forcePairs
  unsigned int host_npart = 0u;

  // pairs==-1: un-initialized device memory
  // pairs==0: return forces using atomicAdd
  // pairs==2: return forces using a global memory reduction
  int pairs = -1;
  // variables for forces and energies calculated pre-reduction
  float *dev_forcePairs = nullptr;
  float *dev_energyBlocks = nullptr;

  // run a single force calculation and return the time it takes using
  // high-precision CUDA timers
  float force_benchmark();

  void setup();
  void modpsi_init();
  void set_params(double boxz, double prefactor, double maxPWerror,
                  double far_switch_radius, int bessel_cutoff);
  void tune(double maxPWerror, double far_switch_radius, int bessel_cutoff);
  void sanity_checks_periodicity() const;
  void sanity_checks_cell_structure() const;
};

#endif // MMM1D_GPU
#endif
