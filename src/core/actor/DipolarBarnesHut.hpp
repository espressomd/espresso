/*
 * Copyright (C) 2016-2019 The ESPResSo project
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
#ifndef ACTOR_DIPOLARBARNESHUT_HPP
#define ACTOR_DIPOLARBARNESHUT_HPP

#include "config.hpp"

#ifdef DIPOLAR_BARNES_HUT

#include "Actor.hpp"
#include "DipolarBarnesHut_cuda.cuh"
#include "SystemInterface.hpp"

class DipolarBarnesHut : public Actor {
public:
  DipolarBarnesHut(SystemInterface &s);
  ~DipolarBarnesHut() override { deallocBH(&m_bh_data); }

  void set_params(float epssq, float itolsq);

  void computeForces(SystemInterface &s) override;
  void computeEnergy(SystemInterface &s) override;

  void activate();
  void deactivate();

private:
  float m_prefactor;
  float m_epssq;
  float m_itolsq;
  /// Container for pointers to device memory.
  BHData m_bh_data = {0,       0,       0,       nullptr, nullptr,
                      nullptr, nullptr, nullptr, nullptr, nullptr,
                      nullptr, nullptr, nullptr, nullptr};
  int initialize_data_structure(SystemInterface &s);
};

#endif // DIPOLAR_BARNES_HUT

#endif // ACTOR_DIPOLARBARNESHUT_HPP
