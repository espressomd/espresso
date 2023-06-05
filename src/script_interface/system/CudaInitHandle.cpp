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

#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace ScriptInterface {
namespace System {

CudaInitHandle::CudaInitHandle() {
  add_parameters({
#ifdef CUDA
      {"device",
       [this](Variant const &v) {
         if (context()->is_head_node()) {
           cuda_set_device(get_value<int>(v));
         }
       },
       [this]() {
         return (context()->is_head_node()) ? cuda_get_device() : 0;
       }},
#endif // CUDA
  });
}

#ifdef CUDA
/**
 * @brief Silently ignore CUDA exceptions.
 * This is useful when querying the properties of CUDA devices that may
 * not have a suitable CUDA version, or when there is no compatible CUDA
 * device available.
 */
template <typename F> static void skip_cuda_errors(F &&fun) {
  try {
    fun();
  } catch (cuda_runtime_error const &) {
  }
}
#endif // CUDA

Variant CudaInitHandle::do_call_method(std::string const &name,
                                       VariantMap const &parameters) {
  if (name == "list_devices") {
    std::unordered_map<int, std::string> devices{};
#ifdef CUDA
    if (context()->is_head_node()) {
      // only GPUs on the head node can be used
      auto n_gpus = 0;
      skip_cuda_errors([&n_gpus]() { n_gpus = cuda_get_n_gpus(); });
      for (int i = 0; i < n_gpus; ++i) {
        skip_cuda_errors([&devices, i]() {
          char gpu_name_buffer[4 + 64];
          cuda_get_gpu_name(i, gpu_name_buffer);
          devices[i] = std::string{gpu_name_buffer};
        });
      }
    }
#endif // CUDA
    return make_unordered_map_of_variants(devices);
  }
  if (name == "list_devices_properties") {
    std::unordered_map<std::string, std::unordered_map<int, Variant>> dict{};
#ifdef CUDA
    std::vector<EspressoGpuDevice> devices = cuda_gather_gpus();
    for (auto const &dev : devices) {
      auto const hostname = std::string{dev.proc_name};
      if (dict.count(hostname) == 0) {
        dict[hostname] = {};
      }
      std::unordered_map<std::string, Variant> dev_properties = {
          {"name", std::string{dev.name}},
          {"compute_capability",
           Variant{std::vector<int>{
               {dev.compute_capability_major, dev.compute_capability_minor}}}},
          {"cores", dev.n_cores},
          {"total_memory", dev.total_memory},
      };
      dict[hostname][dev.id] = std::move(dev_properties);
    }
#endif // CUDA
    return make_unordered_map_of_variants(dict);
  }
  if (name == "get_n_gpus") {
    auto n_gpus = 0;
#ifdef CUDA
    if (context()->is_head_node()) {
      // only GPUs on the head node can be used
      skip_cuda_errors([&n_gpus]() { n_gpus = cuda_get_n_gpus(); });
    }
#endif // CUDA
    return n_gpus;
  }
  return {};
}

} // namespace System
} // namespace ScriptInterface
