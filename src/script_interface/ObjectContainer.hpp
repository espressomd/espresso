/*
 * Copyright (C) 2023 The ESPResSo project
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
#ifndef SCRIPT_INTERFACE_OBJECT_CONTAINER_HPP
#define SCRIPT_INTERFACE_OBJECT_CONTAINER_HPP

#include "script_interface/auto_parameters/AutoParameters.hpp"

#include <type_traits>

namespace ScriptInterface {

/**
 * @brief Base class for containers whose @c BaseType might be a full
 * specialization of @ref AutoParameters.
 */
template <template <typename...> class Container, typename ManagedType,
          class BaseType,
          class =
              std::enable_if_t<std::is_base_of_v<ObjectHandle, ManagedType>>>
using ObjectContainer = std::conditional_t<
    std::is_same_v<BaseType, ObjectHandle>,
    AutoParameters<Container<ManagedType, BaseType>, BaseType>, BaseType>;

} // namespace ScriptInterface

#endif
