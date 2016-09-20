/*
  Copyright (C) 2010,2012,2013,2014,2015,2016 The ESPResSo project
  
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

#include "config.hpp"
#include "utils.hpp"
#include "cuda_init.hpp"
#include "communication.hpp"

#include <mpi.h>
#include <string.h>
#include <set>
#include <iterator>

#ifdef CUDA

/** Helper class force device set 
 */

struct CompareDevices {
  bool operator()(const EspressoGpuDevice &a, const EspressoGpuDevice &b) const {
    const int name_comp = strncmp(a.proc_name, b.proc_name, 63);
    /* Both devs are from the same node, order by id */
    if(name_comp == 0)
      return a.id < b.id;
    else
      return name_comp < 0;
  }
};

/** Gather list of CUDA devices on all nodes on the master node
    It relies on MPI_Get_processor_name() to get a unique identifier of
    the physical node, as opposed to the logical rank of which there can
    be more than one on one node.
 */

std::vector<EspressoGpuDevice> cuda_gather_gpus(void) {
  int n_gpus = cuda_get_n_gpus();
  char proc_name[MPI_MAX_PROCESSOR_NAME];
  int proc_name_len;
  /* List of local devices */
  std::vector<EspressoGpuDevice> devices;
  /* Global unique device list (only relevant on master) */
  std::vector<EspressoGpuDevice> g_devices;
  int *n_gpu_array = 0;

  MPI_Get_processor_name(proc_name, &proc_name_len);

  /* Truncate to 63 chars to fit struct. */
  if(strlen(proc_name) > 63)
    proc_name[63] = 0;

  for(int i = 0; i < n_gpus; ++i) {
    /* Check if device has at least mininum compute capability */
    if(cuda_check_gpu(i) == ES_OK) {
      EspressoGpuDevice device;
      if(cuda_get_device_props(i, device) == ES_OK){
	strncpy(device.proc_name, proc_name, 64);
	devices.push_back(device);
      }
    }
  }
  
  /** Update n_gpus to number of usable devices */
  n_gpus = devices.size();

  if(this_node == 0) {
    std::set<EspressoGpuDevice, CompareDevices> device_set;
    n_gpu_array = new int[n_nodes];
    MPI_Gather(&n_gpus, 1, MPI_INT, n_gpu_array, 1, MPI_INT, 0, MPI_COMM_WORLD);

    /* insert local devices */
    std::copy(devices.begin(), devices.end(), std::inserter(device_set, device_set.begin()));

    EspressoGpuDevice device;      
    MPI_Status s;
    /* Get devices from other nodes */
    for(int i = 1; i < n_nodes; ++i) {
      for(int j = 0; j < n_gpu_array[i]; ++j) {
	MPI_Recv(&device, sizeof(EspressoGpuDevice), MPI_BYTE, i, 0, MPI_COMM_WORLD, &s);
	device_set.insert(device);
      }      
    }
    /* Copy unique devices to result, if any */
    std::copy(device_set.begin(), device_set.end(), std::inserter(g_devices, g_devices.begin()));
    delete n_gpu_array;
  } else {
    /* Send number of devices to master */
    MPI_Gather(&n_gpus, 1, MPI_INT, n_gpu_array, 1, MPI_INT, 0, MPI_COMM_WORLD);
    /* Send devices to maser */
    for(std::vector<EspressoGpuDevice>::iterator device = devices.begin();
	device != devices.end(); ++device) {
      MPI_Send(&(*device), sizeof(EspressoGpuDevice), MPI_BYTE, 0, 0, MPI_COMM_WORLD);
    }
  }  
  return g_devices;
}

#endif /* CUDA */

