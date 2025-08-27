/*
 * Copyright (C) 2013-2022 The ESPResSo project
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

#include "CudaInitHandle.hpp"

#include "config/config.hpp"

#include "core/cuda/init.hpp"
#include "core/cuda/utils.hpp"

#if defined(ESPRESSO_CUDA) && defined(ESPRESSO_WALBERLA)
#include "walberla_bridge/lattice_boltzmann/lb_walberla_init.hpp"
#endif

#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace ScriptInterface {
namespace System {

CudaInitHandle::CudaInitHandle() {
  add_parameters({
#ifdef ESPRESSO_CUDA
      {"device",
       [this](Variant const &v) {
         if (context()->is_head_node()) {
           cuda_set_device(get_value<int>(v));
         }
       },
       [this]() {
         return (context()->is_head_node()) ? cuda_get_device() : 0;
       }},
#endif // ESPRESSO_CUDA
  });
}

Variant CudaInitHandle::do_call_method(std::string const &name,
                                       VariantMap const &parameters) {
  if (name == "list_devices") {
    std::unordered_map<int, std::string> devices{};
#ifdef ESPRESSO_CUDA
    if (context()->is_head_node()) {
      // only GPUs on the head node can be displayed
      auto n_gpus = 0;
      invoke_skip_cuda_exceptions([&n_gpus]() { n_gpus = cuda_get_n_gpus(); });
      for (int i = 0; i < n_gpus; ++i) {
        invoke_skip_cuda_exceptions(
            [&devices, i]() { devices[i] = cuda_get_gpu_name(i); });
      }
    }
#endif // ESPRESSO_CUDA
    return make_unordered_map_of_variants(devices);
  }
  if (name == "list_devices_properties") {
    std::unordered_map<std::string, std::unordered_map<int, Variant>> dict{};
#ifdef ESPRESSO_CUDA
    std::vector<EspressoGpuDevice> devices = cuda_gather_gpus();
    for (auto const &dev : devices) {
      auto const hostname = dev.proc_name;
      if (not dict.contains(hostname)) {
        dict[hostname] = {};
      }
      std::unordered_map<std::string, Variant> dev_properties = {
          {"name", dev.name},
          {"compute_capability",
           Variant{std::vector<int>{
               {dev.compute_capability_major, dev.compute_capability_minor}}}},
          {"cores", dev.n_cores},
          {"total_memory", dev.total_memory},
      };
      dict[hostname][dev.id] = std::move(dev_properties);
    }
#endif // ESPRESSO_CUDA
    return make_unordered_map_of_variants(dict);
  }
  if (name == "get_n_gpus") {
    auto n_gpus = 0;
#ifdef ESPRESSO_CUDA
    auto const devices = cuda_gather_gpus();
    n_gpus = static_cast<int>(devices.size());
#endif // ESPRESSO_CUDA
    return n_gpus;
  }
#if defined(ESPRESSO_CUDA) && defined(ESPRESSO_WALBERLA)
  if (name == "set_device_id_per_rank") {
    if (cuda_get_n_gpus()) {
      set_device_id_per_rank();
    }
    return {};
  }
#endif
  return {};
}

} // namespace System
} // namespace ScriptInterface
