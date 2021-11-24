/*
 * Copyright (C) 2021 The ESPResSo project
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
#ifndef ESPRESSO_OBJECT_CONTAINER_MPI_GUARD_HPP
#define ESPRESSO_OBJECT_CONTAINER_MPI_GUARD_HPP

#include <boost/utility/string_ref.hpp>

#include <cstddef>

/**
 * @brief Prevent object container serialization.
 *
 * The @ref ScriptInterface::ObjectHandle framework doesn't support
 * recursive deserialization. When an object container such as
 * @ref ScriptInterface::ObjectList is deserialized, the contained
 * objects are deserialized on the head node only, which leads to
 * silent bugs in simulations.
 *
 * This function needs to be called from an object container
 * <tt>get_internal_state()</tt> method to throw a runtime error
 * when the container is not empty and the MPI world size is
 * greater than 1.
 *
 * @param name        Name of the object container
 * @param n_elements  Number of elements in the container
 */
void object_container_mpi_guard(boost::string_ref const &name,
                                std::size_t n_elements);

#endif
