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

#pragma once

#include <list>
#include <memory>
#include <stack>

/**
 * @brief Manager to control the lifetime of shared resources.
 *
 * Resources that need to be available globally, for example
 * via singletons, need to be destructed after all objects
 * that depend on them have already been destructed.
 * It is not safe to rely on static global variables, since the
 * order of destruction is not guaranteed across translation units.
 * When such static objects reside in different translation units,
 * they can expire in any order, potentially creating race conditions
 * if one static object relies on the other it its destructor.
 *
 * This class "locks" resources by storing a shared pointer to them,
 * ensuring that the resources lifetime is extended by the lifetime
 * of the class instance. Consumer code can then keep this class
 * instance alive until the resources are no longer needed, at which
 * point the class instance can be released. Type erasure is used to
 * hide the implementation details of the resources being locked.
 *
 * Multiple resources can be locked, using a LIFO (last-in, first-out)
 * container to ensure that the resources are freed in a controlled order.
 * This behavior is critical to avoid undefined behavior when one global
 * resource's destruction depends on the existence of another global
 * resource, and cannot be achieved with STL containers like @p std::stack
 * or @p std::vector, for which the destruction order of the stored data
 * is under-specified.
 */
class ResourceManager {
  class ResourceLock {
  public:
    virtual ~ResourceLock() = default;
  };

  template <typename T> class ResourceLockImpl : public ResourceLock {
    std::shared_ptr<T> m_resource;

  public:
    explicit ResourceLockImpl(std::shared_ptr<T> const &resource)
        : m_resource(resource) {}
    ~ResourceLockImpl() { m_resource.reset(); }
  };

  template <typename T> using LifoList = std::stack<T, std::list<T>>;

  LifoList<std::unique_ptr<ResourceLock>> m_resources;

public:
  ResourceManager() = default;

  ~ResourceManager() {
    while (not m_resources.empty()) {
      m_resources.pop();
    }
  }

  template <typename T> void acquire_lock(std::shared_ptr<T> resource) {
    m_resources.emplace(std::make_unique<ResourceLockImpl<T>>(resource));
  }
};
