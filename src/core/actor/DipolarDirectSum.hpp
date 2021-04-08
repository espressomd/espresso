/*
 * Copyright (C) 2010-2019 The ESPResSo project
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
#ifndef ACTOR_DIPOLARDIRECTSUM_HPP
#define ACTOR_DIPOLARDIRECTSUM_HPP

#include "config.hpp"

#ifdef DIPOLAR_DIRECT_SUM

#include "Actor.hpp"
#include "SystemInterface.hpp"
#include "cuda_interface.hpp"
#include "electrostatics_magnetostatics/dipole.hpp"
#include "errorhandling.hpp"
#include "grid.hpp"

#include <memory>

void DipolarDirectSum_kernel_wrapper_energy(float k, int n, float *pos,
                                            float *dip, float box_l[3],
                                            int periodic[3], float *E);
void DipolarDirectSum_kernel_wrapper_force(float k, int n, float *pos,
                                           float *dip, float *f, float *torque,
                                           float box_l[3], int periodic[3]);

class DipolarDirectSum : public Actor {
public:
  DipolarDirectSum(SystemInterface &s) {

    k = static_cast<float>(dipole.prefactor);

    if (!s.requestFGpu())
      runtimeErrorMsg() << "DipolarDirectSum needs access to forces on GPU!";

    if (!s.requestRGpu())
      runtimeErrorMsg() << "DipolarDirectSum needs access to positions on GPU!";

    if (!s.requestDipGpu())
      runtimeErrorMsg() << "DipolarDirectSum needs access to dipoles on GPU!";
  };
  void computeForces(SystemInterface &s) override {
    float box[3];
    int per[3];
    for (int i = 0; i < 3; i++) {
      box[i] = static_cast<float>(s.box()[i]);
      per[i] = (box_geo.periodic(i));
    }
    DipolarDirectSum_kernel_wrapper_force(
        k, static_cast<int>(s.npart_gpu()), s.rGpuBegin(), s.dipGpuBegin(),
        s.fGpuBegin(), s.torqueGpuBegin(), box, per);
  };
  void computeEnergy(SystemInterface &s) override {
    float box[3];
    int per[3];
    for (int i = 0; i < 3; i++) {
      box[i] = static_cast<float>(s.box()[i]);
      per[i] = (box_geo.periodic(i));
    }
    DipolarDirectSum_kernel_wrapper_energy(
        k, static_cast<int>(s.npart_gpu()), s.rGpuBegin(), s.dipGpuBegin(), box,
        per, (&(((CUDA_energy *)s.eGpu())->dipolar)));
  };

private:
  float k;
};

void activate_dipolar_direct_sum_gpu();
void deactivate_dipolar_direct_sum_gpu();

extern std::unique_ptr<DipolarDirectSum> dipolarDirectSum;

#endif

#endif
