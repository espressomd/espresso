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

/**
 * Public interface for the ESPResSo system and device memory management.
 *
 * This interface is responsible for device memory management and memory
 * layout conversion (AoS to SoA conversion). On the CPU, particle data
 * is stored in an Array of Structures (AoS). On the GPU, some algorithms
 * have better performance with a Structure of Arrays (SoA). When data
 * is synchronized between host and device memory, particle objects are
 * copied to device memory in the form of an AoS. A few particle properties
 * are further copied ("split") into a SoA by the split callback.
 *
 * This interface provides getters in the form <tt>xGpuBegin()</tt> and
 * <tt>xGpuEnd()</tt> (with @c x a particle property, e.g. @c r for position)
 * that return pointers to device memory for each array in the SoA.
 * To register a new particle property in the split callback, call the
 * corresponding <tt>requestXGpu()</tt> method, with @c X the uppercase
 * version of the particle property @c x. It is currently impossible to
 * undo this action.
 *
 * Since properties in the @ref Particle struct depend on the config file,
 * it is necessary to provide a means to check if it is even possible to
 * carry out the split. This is achieved by the <tt>hasXGpu()</tt> methods.
 * The same value is returned from <tt>requestXGpu()</tt> and it needs to
 * be checked at the call site to throw an error. This class won't throw
 * exceptions.
 */
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
