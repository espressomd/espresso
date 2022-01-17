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
#ifndef _CUDA_INIT_H
#define _CUDA_INIT_H

#include "config.hpp"

#ifdef CUDA

#include <cstddef>
#include <vector>

/** Struct to hold information relevant to ESPResSo
 *  about GPUs. Should contain only fixed-length plain
 *  old datatypes, as it is intended for MPI communication.
 */
struct EspressoGpuDevice {
  /** Local CUDA device id */
  int id;
  /** Local CUDA device name */
  char name[64];
  /** Node identification */
  char proc_name[64];
  /** MPI process identification */
  int node;
  /** Compute capability (major) */
  int compute_capability_major;
  /** Compute capability (minor) */
  int compute_capability_minor;
  /** Total Memory */
  std::size_t total_memory;
  /** Number of cores */
  int n_cores;
};

/** Initializes the CUDA stream.
 */
void cuda_init();

/** Get the number of CUDA devices.
 *
 *  @return the number of GPUs.
 */
int cuda_get_n_gpus();

/** Check that a given GPU has compute capability.
 *  The minimal compute capability required by ESPResSo is
 *  \ref computeCapabilityMinMajor . \ref computeCapabilityMinMinor .
 *
 *  @param dev CUDA device number
 *  @return \ref ES_OK if the GPU meets the requirements, else \ref ES_ERROR.
 */
int cuda_check_gpu_compute_capability(int dev);

/** Get the name of a CUDA device.
 *
 *  @param[in]  dev the CUDA device number to ask the name for
 *  @param[out] name a buffer to write the name to, at least 64 characters
 */
void cuda_get_gpu_name(int dev, char name[64]);

/** Choose a device for future CUDA computations.
 *
 *  @param dev the device to use
 */
void cuda_set_device(int dev);

/** Get the current CUDA device.
 *
 *  @return the current device's number.
 */
int cuda_get_device();

/** Test if actual CUDA device works.
 *  @return \ref ES_OK on success, \ref ES_ERROR else.
 */
int cuda_test_device_access();

/**
 * Check that a device is available, that its compute capability
 * is sufficient for ESPResSo, and that data can be written to
 * and read from it. Otherwise, throw an exception.
 */
void cuda_check_device();

/** Gather unique list of CUDA devices on all nodes.
 *  @return vector of device properties.
 */
std::vector<EspressoGpuDevice> cuda_gather_gpus();

/** Get properties of a CUDA device
 *  @param dev CUDA device number
 */
EspressoGpuDevice cuda_get_device_props(int dev);

#endif // ifdef CUDA
#endif
