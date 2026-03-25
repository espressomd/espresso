/*
 * Copyright (C) 2010-2026 The ESPResSo project
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

#include "config/config.hpp"

#ifdef ESPRESSO_CUDA

#include "init.hpp"
#include "utils.hpp"

#include "communication.hpp"

#include <utils/mpi/gather_buffer.hpp>

#include <mpi.h>

#include <algorithm>
#include <iterator>
#include <set>
#include <string>
#include <utility>
#include <vector>

/** Helper class for device sets.
 */
struct CompareDevices {
  bool operator()(const EspressoGpuDevice &a,
                  const EspressoGpuDevice &b) const {
    auto const name_comp = a.proc_name.compare(b.proc_name);
    /* if both devices are from the same node, order by id */
    return (name_comp == 0) ? a.id < b.id : name_comp < 0;
  }
};

/** Gather list of CUDA devices from all nodes on the head node.
 *  It relies on <tt>MPI_Get_processor_name()</tt> to get a unique identifier
 *  of the physical node, as opposed to the logical rank of which there can
 *  be more than one per node.
 */
std::vector<EspressoGpuDevice> cuda_gather_gpus() {
  std::vector<EspressoGpuDevice> devices;

  int n_devices = 0;
  invoke_skip_cuda_exceptions(
      [&n_devices]() { n_devices = cuda_get_n_gpus(); });

  invoke_skip_cuda_exceptions([&devices, n_devices]() {
    for (int i = 0; i < n_devices; ++i) {
      devices.emplace_back(cuda_get_device_props(i));
    }
  });

  Utils::Mpi::gather_buffer(devices, ::comm_cart);

  if (this_node == 0) {
    std::set<EspressoGpuDevice, CompareDevices> device_set;
    std::ranges::copy(devices, std::inserter(device_set, device_set.begin()));
    devices.clear();
    std::ranges::copy(device_set, std::back_inserter(devices));
  } else {
    devices.clear();
  }
  return devices;
}

namespace detail {
std::pair<int, std::string> get_node_info() {
  int proc_name_len;
  char proc_name[MPI_MAX_PROCESSOR_NAME];
  MPI_Get_processor_name(proc_name, &proc_name_len);
  return {this_node, proc_name};
}
} // namespace detail

void cuda_on_program_start() {
  if (::communicator.this_node == 0) {
    invoke_skip_cuda_exceptions(cuda_init);
  }
}

#endif // ESPRESSO_CUDA
