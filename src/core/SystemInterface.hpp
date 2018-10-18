/*
  Copyright (C) 2014-2018 The ESPResSo project

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
#ifndef SYSTEMINTERFACE_H
#define SYSTEMINTERFACE_H

#include "Vector.hpp"
#include "config.hpp"
#include <vector>

/** @todo: Turn needsXY in getter/setter **/

class SystemInterface {
public:
  SystemInterface()
      : m_needsRGpu(false), m_needsVGpu(false), m_needsQGpu(false),
        m_needsDirectorGpu(false), m_needsFGpu(false), m_needsDipGpu(false),
        m_needsTorqueGpu(false){};
  typedef Vector3d Vector3;
  typedef double Real;

  virtual void init() = 0;
  virtual void update() = 0;

  virtual float *rGpuBegin() { return 0; };
  virtual float *rGpuEnd() { return 0; };
  virtual bool hasRGpu() { return false; };
  virtual bool requestRGpu() {
    m_needsRGpu = hasRGpu();
    return m_needsRGpu;
  }

  virtual float *dipGpuBegin() { return 0; };
  virtual float *dipGpuEnd() { return 0; };
  virtual bool hasDipGpu() { return false; };
  virtual bool requestDipGpu() {
    m_needsDipGpu = hasDipGpu();
    return m_needsDipGpu;
  }

  virtual float *torqueGpuBegin() { return 0; };
  virtual float *torqueGpuEnd() { return 0; };
  virtual bool hasTorqueGpu() { return false; };
  virtual bool requestTorqueGpu() {
    m_needsTorqueGpu = hasTorqueGpu();
    return m_needsTorqueGpu;
  }

  virtual float *vGpuBegin() { return 0; };
  virtual float *vGpuEnd() { return 0; };
  virtual bool hasVGpu() { return false; };
  virtual bool requestVGpu() {
    m_needsVGpu = hasVGpu();
    return m_needsVGpu;
  }

  virtual float *qGpuBegin() { return 0; };
  virtual float *qGpuEnd() { return 0; };
  virtual bool hasQGpu() { return false; };
  virtual bool requestQGpu() {
    m_needsQGpu = hasQGpu();
    return m_needsQGpu;
  }

  virtual float *fGpuBegin() { return 0; };
  virtual float *fGpuEnd() { return 0; };
  virtual float *eGpu() { return 0; };
  virtual bool hasFGpu() { return false; };
  virtual bool requestFGpu() {
    m_needsFGpu = hasFGpu();
    return m_needsFGpu;
  }

  virtual float *directorGpuBegin() { return 0; };
  virtual float *directorGpuEnd() { return 0; };
  virtual bool hasDirectorGpu() { return false; };
  virtual bool requestDirectorGpu() {
    m_needsDirectorGpu = hasDirectorGpu();
    return m_needsDirectorGpu;
  }

  virtual unsigned int npart_gpu() { return 0; };
  virtual Vector3 box() const = 0;

  virtual bool needsRGpu() { return m_needsRGpu; };
  virtual bool needsDipGpu() { return m_needsRGpu; };
  virtual bool needsQGpu() { return m_needsQGpu; };
  virtual bool needsDirectorGpu() { return m_needsDirectorGpu; };
  virtual bool needsFGpu() { return m_needsFGpu; };
  virtual bool needsTorqueGpu() { return m_needsTorqueGpu; };
  virtual ~SystemInterface() = default;

protected:
  bool m_needsRGpu;
  bool m_needsVGpu;
  bool m_needsQGpu;
  bool m_needsDirectorGpu;
  bool m_needsFGpu;
  bool m_needsDipGpu;
  bool m_needsTorqueGpu;
};

#endif
