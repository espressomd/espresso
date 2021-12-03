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
#ifndef SYSTEMINTERFACE_H
#define SYSTEMINTERFACE_H

#include "config.hpp"

#include <utils/Vector.hpp>

#include <cstddef>

class SystemInterface {
public:
  SystemInterface() = default;
  virtual ~SystemInterface() = default;

  virtual void init() = 0;
  virtual void update() = 0;

  virtual float *rGpuBegin() { return nullptr; };
  virtual float *rGpuEnd() { return nullptr; };
  virtual bool hasRGpu() { return false; };
  virtual bool requestRGpu() { return false; }

  virtual float *dipGpuBegin() { return nullptr; };
  virtual float *dipGpuEnd() { return nullptr; };
  virtual bool hasDipGpu() { return false; };
  virtual bool requestDipGpu() { return false; }

  virtual float *torqueGpuBegin() { return nullptr; };
  virtual float *torqueGpuEnd() { return nullptr; };
  virtual bool hasTorqueGpu() { return false; };
  virtual bool requestTorqueGpu() { return false; }

  virtual float *qGpuBegin() { return nullptr; };
  virtual float *qGpuEnd() { return nullptr; };
  virtual bool hasQGpu() { return false; };
  virtual bool requestQGpu() { return false; }

  virtual float *fGpuBegin() { return nullptr; };
  virtual float *fGpuEnd() { return nullptr; };
  virtual bool hasFGpu() { return false; };
  virtual bool requestFGpu() { return false; }

  virtual float *eGpu() { return nullptr; };

  virtual std::size_t npart_gpu() const { return 0; };
  virtual Utils::Vector3d box() const = 0;
};

#endif
