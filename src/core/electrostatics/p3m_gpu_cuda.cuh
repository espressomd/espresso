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

#pragma once

#include "system/GpuParticleData.hpp"

#include <utils/Vector.hpp>

#include <cstddef>
#include <memory>

struct P3MGpuParams;

void p3m_gpu_init(std::shared_ptr<P3MGpuParams> &p3m_gpu_data_ptr, int cao,
                  Utils::Vector3i const &mesh, double alpha,
                  Utils::Vector3d const &box_l, std::size_t n_part);
void p3m_gpu_add_farfield_force(P3MGpuParams &data, GpuParticleData &gpu,
                                double prefactor, std::size_t n_part);
