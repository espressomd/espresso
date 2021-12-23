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
#ifndef ESPRESSOSYSTEMINTERFACE_H
#define ESPRESSOSYSTEMINTERFACE_H

#include "SystemInterface.hpp"
#include "config.hpp"
#include "cuda_interface.hpp"

#include <utils/Vector.hpp>

#include <cstddef>

/**
 * @brief CUDA implementation of @ref SystemInterface.
 *
 * When data is synchronized between host and device memory, a subset
 * of the @ref Particle struct is copied from each particle on the host
 * to the corresponding @ref CUDA_particle_data struct on the device via
 * @ref EspressoSystemInterface::gatherParticles(). Once the transfer is
 * complete, the particle AoS on the device is copied (or "split") to
 * a SoA via @ref EspressoSystemInterface::split_particle_struct().
 */
class EspressoSystemInterface : public SystemInterface {
public:
  EspressoSystemInterface() = default;
  ~EspressoSystemInterface() override = default;

  static EspressoSystemInterface &Instance() {
    if (!m_instance)
      m_instance = new EspressoSystemInterface;

    return *m_instance;
  };

  void init() override;
  void update() override;

#ifdef CUDA
  float *rGpuBegin() override { return m_r_gpu_begin; };
  bool hasRGpu() override { return true; };
  void requestRGpu() override {
    m_needsRGpu = hasRGpu();
    m_splitParticleStructGpu |= m_needsRGpu;
    m_gpu |= m_needsRGpu;
    enableParticleCommunication();
  };

#ifdef DIPOLES
  float *dipGpuBegin() override { return m_dip_gpu_begin; };
  bool hasDipGpu() override { return true; };
  void requestDipGpu() override {
    m_needsDipGpu = hasDipGpu();
    m_splitParticleStructGpu |= m_needsRGpu;
    m_gpu |= m_needsRGpu;
    enableParticleCommunication();
  };
#endif

#ifdef ELECTROSTATICS
  float *qGpuBegin() override { return m_q_gpu_begin; };
  bool hasQGpu() override { return true; };
  void requestQGpu() override {
    m_needsQGpu = hasQGpu();
    m_splitParticleStructGpu |= m_needsQGpu;
    m_gpu |= m_needsQGpu;
    enableParticleCommunication();
  };
#endif

  void requestParticleStructGpu() {
    m_needsParticleStructGpu = true;
    m_gpu |= m_needsParticleStructGpu;
    enableParticleCommunication();
  };

  float *fGpuBegin() override { return gpu_get_particle_force_pointer(); };
  bool hasFGpu() override { return true; };
  void requestFGpu() override {
    m_needsFGpu = hasFGpu();
    m_gpu |= m_needsFGpu;
    enableParticleCommunication();
  };

#ifdef ROTATION
  float *torqueGpuBegin() override {
    return gpu_get_particle_torque_pointer();
  };
  bool hasTorqueGpu() override { return true; };
  void requestTorqueGpu() override {
    m_needsTorqueGpu = hasTorqueGpu();
    m_gpu |= m_needsTorqueGpu;
    enableParticleCommunication();
  };
#endif

  float *eGpu() override {
    // cast pointer from struct of floats to array of floats
    // https://stackoverflow.com/a/29278260
    return reinterpret_cast<float *>(gpu_get_energy_pointer());
  };

#endif // ifdef CUDA

  Utils::Vector3d box() const override;

  std::size_t npart_gpu() const override {
#ifdef CUDA
    return gpu_get_particle_pointer().size();
#else
    return 0;
#endif
  };

protected:
  static EspressoSystemInterface *m_instance;

  void gatherParticles();
  void split_particle_struct();
#ifdef CUDA
  void enableParticleCommunication();
  void reallocDeviceMemory(std::size_t n);

private:
  std::size_t m_gpu_npart = 0;
  bool m_gpu = false;

  float *m_r_gpu_begin = nullptr;
  float *m_dip_gpu_begin = nullptr;
  float *m_q_gpu_begin = nullptr;

  bool m_needsParticleStructGpu = false;
  bool m_splitParticleStructGpu = false;

  bool m_needsRGpu = false;
  bool m_needsQGpu = false;
  bool m_needsFGpu = false;
  bool m_needsDipGpu = false;
  bool m_needsTorqueGpu = false;
#endif
};

#endif
