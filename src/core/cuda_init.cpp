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

#include "config.hpp"

#ifdef CUDA

#include "cuda_init.hpp"
#include "cuda_utils.hpp"

#include "communication.hpp"
#include "errorhandling.hpp"

#include <utils/constants.hpp>

#include <mpi.h>

#include <algorithm>
#include <cstring>
#include <iterator>
#include <set>
#include <vector>

/** Helper class for device sets.
 */
struct CompareDevices {
  bool operator()(const EspressoGpuDevice &a,
                  const EspressoGpuDevice &b) const {
    const int name_comp = strncmp(a.proc_name, b.proc_name, 63);
    /* if both devices are from the same node, order by id */
    if (name_comp == 0)
      return a.id < b.id;

    return name_comp < 0;
  }
};

/** Gather list of CUDA devices on all nodes on the head node.
 *  It relies on <tt>MPI_Get_processor_name()</tt> to get a unique identifier
 *  of the physical node, as opposed to the logical rank of which there can
 *  be more than one per node.
 */
static std::vector<EspressoGpuDevice> mpi_cuda_gather_gpus_local() {
  /* List of local devices */
  std::vector<EspressoGpuDevice> devices_local;
  /* Global unique device list (only relevant on the head node) */
  std::vector<EspressoGpuDevice> devices_global;

  int n_devices;
  try {
    n_devices = cuda_get_n_gpus();
  } catch (cuda_runtime_error const &err) {
    n_devices = 0;
  }

  int proc_name_len;
  char proc_name[MPI_MAX_PROCESSOR_NAME];
  MPI_Get_processor_name(proc_name, &proc_name_len);
  if (std::strlen(proc_name) > 63)
    proc_name[63] = '\0';

  for (int i = 0; i < n_devices; ++i) {
    try {
      EspressoGpuDevice device = cuda_get_device_props(i);
      std::strncpy(device.proc_name, proc_name, 64);
      device.proc_name[63] = '\0';
      device.node = this_node;
      devices_local.push_back(device);
    } catch (cuda_runtime_error const &err) {
      // pass
    }
  }

  int const n_gpus = static_cast<int>(devices_local.size());

  if (this_node == 0) {
    std::set<EspressoGpuDevice, CompareDevices> device_set;
    int *n_gpu_array = new int[n_nodes];
    MPI_Gather(&n_gpus, 1, MPI_INT, n_gpu_array, 1, MPI_INT, 0, MPI_COMM_WORLD);

    /* insert local devices */
    std::copy(devices_local.begin(), devices_local.end(),
              std::inserter(device_set, device_set.begin()));

    EspressoGpuDevice device;
    MPI_Status s;
    /* Get devices from other nodes */
    for (int i = 1; i < n_nodes; ++i) {
      for (int j = 0; j < n_gpu_array[i]; ++j) {
        MPI_Recv(&device, sizeof(EspressoGpuDevice), MPI_BYTE, i, 0,
                 MPI_COMM_WORLD, &s);
        device_set.insert(device);
      }
    }
    /* Copy unique devices to result, if any */
    std::copy(device_set.begin(), device_set.end(),
              std::inserter(devices_global, devices_global.begin()));
    delete[] n_gpu_array;
  } else {
    /* Send number of devices to head node */
    MPI_Gather(&n_gpus, 1, MPI_INT, nullptr, 1, MPI_INT, 0, MPI_COMM_WORLD);
    /* Send devices to head node */
    for (auto const &device : devices_local) {
      MPI_Send(&device, sizeof(EspressoGpuDevice), MPI_BYTE, 0, 0,
               MPI_COMM_WORLD);
    }
  }
  return devices_global;
}

REGISTER_CALLBACK_MAIN_RANK(mpi_cuda_gather_gpus_local)

std::vector<EspressoGpuDevice> cuda_gather_gpus() {
  return mpi_call(Communication::Result::main_rank, mpi_cuda_gather_gpus_local);
}

#endif /* CUDA */
