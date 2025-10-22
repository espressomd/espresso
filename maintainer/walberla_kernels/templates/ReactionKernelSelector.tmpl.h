/*
 * Copyright (C) 2022-2023 The ESPResSo project
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

#pragma once

{% for i in range(1, max_num_reactants + 1) %}
#include "{{class_name}}_{{i}}_{{precision_suffix[True]}}.h"
#include "{{class_name}}_{{i}}_{{precision_suffix[False]}}.h"
{% endfor %}
#if defined(__CUDACC__)
{% for i in range(1, max_num_reactants + 1) %}
#include "{{class_name}}_{{i}}_{{precision_suffix[True]}}_CUDA.h"
#include "{{class_name}}_{{i}}_{{precision_suffix[False]}}_CUDA.h"
{% endfor %}
#endif

#include <domain_decomposition/BlockDataID.h>

#include <cstddef>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

namespace walberla {
namespace detail {
namespace {{class_name}}Selector {

{% for gpu_suffix in ["", "GPU"] %}
{% set cuda_suffix = "" %}
{% set helper_suffix = "" %}
{% if gpu_suffix == "GPU" %}
{% set cuda_suffix = "_CUDA" %}
{% set helper_suffix = "_gpu" %}
#if defined(__CUDACC__)
{% endif %}

template <typename FloatType = double, std::size_t N = 1> struct KernelTrait{{gpu_suffix}} {
  using {{class_name}}{{gpu_suffix}} =
      {{namespace}}::{{class_name}}_1_{{precision_suffix[True]}}{{cuda_suffix}};
};
{% for i in range(2, max_num_reactants + 1) %}
template <> struct KernelTrait{{gpu_suffix}}<double, {{i}}> {
  using {{class_name}}{{gpu_suffix}} =
      {{namespace}}::{{class_name}}_{{i}}_{{precision_suffix[True]}}{{cuda_suffix}};
};
{% endfor %}
{% for i in range(1, max_num_reactants + 1) %}
template <> struct KernelTrait{{gpu_suffix}}<float, {{i}}> {
  using {{class_name}}{{gpu_suffix}} =
      {{namespace}}::{{class_name}}_{{i}}_{{precision_suffix[False]}}{{cuda_suffix}};
};
{% endfor %}

template <typename FloatType, class Reactant, std::size_t... ints>
auto get_kernel_impl{{helper_suffix}}(const std::vector<std::shared_ptr<Reactant>> &reactants,
                     const double coefficient,
                     {% if class_name == 'ReactionKernelIndexed' -%}
                     const BlockDataID &indexFieldID,
                     {% endif -%}
                     std::index_sequence<ints...> int_seq) {
  auto kernel = std::make_shared<
      typename KernelTrait{{gpu_suffix}}<FloatType, int_seq.size()>::{{class_name}}{{gpu_suffix}}>(
      {% if class_name == 'ReactionKernelIndexed' -%}
      indexFieldID,
      {% endif -%}
      walberla::BlockDataID(
          reactants[ints]->get_species()->get_density_id())...,
      numeric_cast<FloatType>(reactants[ints]->get_order())...,
      numeric_cast<FloatType>(coefficient),
      numeric_cast<FloatType>(reactants[ints]->get_stoech_coeff())...);

  std::function<void(IBlock *)> sweep = [kernel](IBlock * b) { kernel->run(b); };
  return sweep;
}

template <typename FloatType, class Reactant, class... Args>
auto get_kernel_impl{{helper_suffix}}(const std::vector<std::shared_ptr<Reactant>> &reactants,
                     Args... args) {
  switch (reactants.size()) {
{% for i in range(1, max_num_reactants + 1) %}
  case {{i}}:
    return get_kernel_impl{{helper_suffix}}<FloatType>(reactants, args...,
                                      std::make_index_sequence<{{i}}>{});
{% endfor %}
  default:
    throw std::runtime_error("reactions of this size are not implemented!");
  }
}

template <class Reactant, class... Args>
auto get_kernel{{helper_suffix}}(const std::vector<std::shared_ptr<Reactant>> &reactants,
                Args... args) {

  const auto is_double_precision =
      reactants[0]->get_species()->is_double_precision();

  if (is_double_precision) {
    return get_kernel_impl{{helper_suffix}}<double>(reactants, args...);
  }

  return get_kernel_impl{{helper_suffix}}<float>(reactants, args...);
}

{% if gpu_suffix == "GPU" %}
#endif
{% endif %}

{% endfor %}

} // namespace {{class_name}}Selector
} // namespace detail
} // namespace walberla
