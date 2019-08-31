/*
  Copyright (C) 2016,2017,2018 The ESPResSo project

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

#include "config.hpp"

#ifdef DIPOLAR_BARNES_HUT

#include "Actor.hpp"
#include "DipolarBarnesHut_cuda.cuh"
#include "EspressoSystemInterface.hpp"
#include "SystemInterface.hpp"
#include "cuda_interface.hpp"
#include "electrostatics_magnetostatics/dipole.hpp"
#include "errorhandling.hpp"
#include "grid.hpp"

#include <iostream>

#ifndef ACTOR_DIPOLARBARNESHUT_HPP
#define ACTOR_DIPOLARBARNESHUT_HPP

// This needs to be done in the .cu file too
typedef float dds_float;

class DipolarBarnesHut : public Actor {
public:
  DipolarBarnesHut(SystemInterface &s, float epssq, float itolsq) {
    k = dipole.prefactor;
    m_epssq = epssq;
    m_itolsq = itolsq;
    setBHPrecision(&m_epssq, &m_itolsq);
    if (!s.requestFGpu())
      std::cerr << "DipolarBarnesHut needs access to forces on GPU!"
                << std::endl;

    if (!s.requestRGpu())
      std::cerr << "DipolarBarnesHut needs access to positions on GPU!"
                << std::endl;

    if (!s.requestDipGpu())
      std::cerr << "DipolarBarnesHut needs access to dipoles on GPU!"
                << std::endl;

    allocBHmemCopy(s.npart_gpu(), &m_bh_data);
  };

  void computeForces(SystemInterface &s) override {
    fillConstantPointers(s.rGpuBegin(), s.dipGpuBegin(), m_bh_data);
    initBHgpu(m_bh_data.blocks);
    buildBoxBH(m_bh_data.blocks);
    buildTreeBH(m_bh_data.blocks);
    summarizeBH(m_bh_data.blocks);
    sortBH(m_bh_data.blocks);
    if (forceBH(&m_bh_data, k, s.fGpuBegin(), s.torqueGpuBegin())) {
      fprintf(
          stderr,
          "forceBH: some of kernels encounter the algorithm functional error");
      errexit();
    }
  };
  void computeEnergy(SystemInterface &s) override {
    fillConstantPointers(s.rGpuBegin(), s.dipGpuBegin(), m_bh_data);
    initBHgpu(m_bh_data.blocks);
    buildBoxBH(m_bh_data.blocks);
    buildTreeBH(m_bh_data.blocks);
    summarizeBH(m_bh_data.blocks);
    sortBH(m_bh_data.blocks);
    if (energyBH(&m_bh_data, k, (&(((CUDA_energy *)s.eGpu())->dipolar)))) {
      fprintf(
          stderr,
          "energyBH: some of kernels encounter the algorithm functional error");
      errexit();
    }
  };

protected:
  float k;
  float m_epssq;
  float m_itolsq;
  BHData m_bh_data = {0,       0,       0,       nullptr, nullptr,
                      nullptr, nullptr, nullptr, nullptr, nullptr,
                      nullptr, nullptr, nullptr, nullptr};
};

void activate_dipolar_barnes_hut(float epssq, float itolsq);
void deactivate_dipolar_barnes_hut();

extern std::unique_ptr<DipolarBarnesHut> dipolarBarnesHut;

#endif

#endif // DIPOLAR_BARNES_HUT
