/*
 * Copyright (C) 2014-2019 The ESPResSo project
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
#ifndef ESPRESSO_CORE_ACTOR_MMM1DGPUFORCE_HPP
#define ESPRESSO_CORE_ACTOR_MMM1DGPUFORCE_HPP

#include "config.hpp"

#ifdef MMM1D_GPU

#include "Actor.hpp"
#include "SystemInterface.hpp"

class Mmm1dgpuForce : public Actor {
public:
  Mmm1dgpuForce(SystemInterface &s);
  ~Mmm1dgpuForce() override;

  // interface methods
  void computeForces(SystemInterface &s) override;
  void computeEnergy(SystemInterface &s) override;
  // configuration methods
  void setup(SystemInterface &s);
  void tune(SystemInterface &s, float _maxPWerror, float _far_switch_radius,
            int _bessel_cutoff);
  void set_params(float _boxz, float _coulomb_prefactor, float _maxPWerror,
                  float _far_switch_radius, int _bessel_cutoff);
  void activate();
  void deactivate();

private:
  // CUDA parameters
  unsigned int numThreads = 64;
  unsigned int numBlocks(SystemInterface const &s) const;

  // the box length currently set on the GPU
  // Needed to make sure it hasn't been modified after inter coulomb was used.
  float host_boxz = 0;
  // the number of particles we had during the last run. Needed to check if we
  // have to realloc dev_forcePairs
  unsigned int host_npart = 0;
  bool need_tune = true;

  // pairs==0: return forces using atomicAdd
  // pairs==1: return force pairs
  // pairs==2: return forces using a global memory reduction
  int pairs = -1;
  // variables for forces and energies calculated pre-reduction
  float *dev_forcePairs = nullptr;
  float *dev_energyBlocks = nullptr;

  // MMM1D parameters
  float coulomb_prefactor = 0;
  float maxPWerror = -1;
  float far_switch_radius = -1;
  int bessel_cutoff = -1;

  // run a single force calculation and return the time it takes using
  // high-precision CUDA timers
  float force_benchmark(SystemInterface &s);

  void modpsi_init();

  // some functions to move MPI dependencies out of the .cu file
  void sanity_checks();
};

#endif
#endif
