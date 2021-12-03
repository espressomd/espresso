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

/** @todo: Turn needsXY in getter/setter **/

class SystemInterface {
public:
  SystemInterface()
      : m_needsRGpu(false), m_needsQGpu(false), m_needsFGpu(false),
        m_needsDipGpu(false), m_needsTorqueGpu(false){};
  typedef Utils::Vector3d Vector3;
  typedef double Real;

  virtual void init() = 0;
  virtual void update() = 0;

  virtual float *rGpuBegin() { return nullptr; };
  virtual float *rGpuEnd() { return nullptr; };
  virtual bool hasRGpu() { return false; };
  virtual bool requestRGpu() {
    m_needsRGpu = hasRGpu();
    return m_needsRGpu;
  }

  virtual float *dipGpuBegin() { return nullptr; };
  virtual float *dipGpuEnd() { return nullptr; };
  virtual bool hasDipGpu() { return false; };
  virtual bool requestDipGpu() {
    m_needsDipGpu = hasDipGpu();
    return m_needsDipGpu;
  }

  virtual float *torqueGpuBegin() { return nullptr; };
  virtual float *torqueGpuEnd() { return nullptr; };
  virtual bool hasTorqueGpu() { return false; };
  virtual bool requestTorqueGpu() {
    m_needsTorqueGpu = hasTorqueGpu();
    return m_needsTorqueGpu;
  }

  virtual float *qGpuBegin() { return nullptr; };
  virtual float *qGpuEnd() { return nullptr; };
  virtual bool hasQGpu() { return false; };
  virtual bool requestQGpu() {
    m_needsQGpu = hasQGpu();
    return m_needsQGpu;
  }

  virtual float *fGpuBegin() { return nullptr; };
  virtual float *fGpuEnd() { return nullptr; };
  virtual float *eGpu() { return nullptr; };
  virtual bool hasFGpu() { return false; };
  virtual bool requestFGpu() {
    m_needsFGpu = hasFGpu();
    return m_needsFGpu;
  }

  virtual std::size_t npart_gpu() const { return 0; };
  virtual Vector3 box() const = 0;

  virtual bool needsRGpu() { return m_needsRGpu; };
  virtual bool needsDipGpu() { return m_needsRGpu; };
  virtual bool needsQGpu() { return m_needsQGpu; };
  virtual bool needsFGpu() { return m_needsFGpu; };
  virtual bool needsTorqueGpu() { return m_needsTorqueGpu; };
  virtual ~SystemInterface() = default;

protected:
  bool m_needsRGpu;
  bool m_needsQGpu;
  bool m_needsFGpu;
  bool m_needsDipGpu;
  bool m_needsTorqueGpu;
};

#endif
