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
#include <type_traits>

/**
 * @brief Queue to deallocate resources before normal program termination.
 *
 * Some resources can only be deallocated when a static global object is
 * still alive, for example a MPI manager, a GPU context or a file handle.
 * This class can be registered as a callback for normal program termination;
 * any registered resource that didn't expire will be forcefully deallocated.
 *
 * To this end, the "client" class needs to implement cleanup methods that
 * deallocate managed resources, to be called by the class destructor when
 * the instance lifetime expires, or by an "attorney" class at normal
 * program termination, whichever comes first. The attorney class creates
 * a callback that is pushed in a queue that gets processed before normal
 * program termination at the user's discretion, or when the queue itself
 * expires, if the client class instance hasn't expired already.
 *
 * Please note the cleanup methods need to be able to run twice, since
 * the client class destructor will be called eventually, possibly after
 * @c __run_exit_handlers() is called. The attorney-client idiom is used to
 * make the private deallocation methods accessible to the cleanup callbacks.
 */
class ResourceCleanup {
  struct Callback {
    virtual ~Callback() = default;
  };

  std::list<std::unique_ptr<Callback>> m_resources;

public:
  ResourceCleanup() = default;
  ~ResourceCleanup() {
    while (not m_resources.empty()) {
      m_resources.pop_back();
    }
  }

  /**
   * @brief Attorney for a resource deallocator.
   * Usage:
   * @code{.cpp}
   * #include "system/ResourceCleanup.hpp"
   * #include <thrust/device_free.h>
   * #include <thrust/device_malloc.h>
   * #include <memory>
   *
   * class MyClass;
   * static std::shared_ptr<ResourceCleanup> resource_queue = nullptr;
   * static std::shared_ptr<MyClass> resource = nullptr;
   *
   * class MyClass {
   * private:
   *   thrust::device_ptr<float> data_dev;
   *   void free_device_memory() { thrust::device_free(data_dev); }
   *   using Cleanup = ResourceCleanup::Attorney<&MyClass::free_device_memory>;
   *   friend Cleanup;
   *
   * public:
   *   MyClass(int n) { data_dev = thrust::device_malloc<float>(n); }
   *   ~MyClass() { free_device_memory(); }
   *   template <class... Args>
   *   static auto make_shared(Args ...args) {
   *     auto obj = std::make_shared<MyClass>(args...);
   *     // comment the next line of code to trigger the CUDA runtime error
   *     // "device free failed: cudaErrorCudartUnloading: driver shutting down"
   *     ::resource_queue->push<MyClass::Cleanup>(obj);
   *     return obj;
   *   }
   * };
   *
   * int main() {
   *   ::resource_queue = std::make_shared<ResourceCleanup>();
   *   ::resource = MyClass::make_shared(5);
   *   // The CUDA primary context expires before normal program termination.
   *   // Yet thrust vectors require a primary context to deallocate device
   *   // memory. Force deallocation before normal program termination.
   *   ::resource_queue.reset();
   * }
   * @endcode
   * @tparam deallocator  Member function that deallocates managed resources.
   */
  template <auto deallocator> class Attorney {
    struct detail {
      template <class C> struct get_class_from_member_function_pointer;
      template <class C, class Ret, class... Args>
      struct get_class_from_member_function_pointer<Ret (C::*)(Args...)> {
        using type = C;
      };
    };
    using Deallocator = decltype(deallocator);
    static_assert(std::is_member_function_pointer_v<Deallocator>);
    using Container =
        typename detail::template get_class_from_member_function_pointer<
            Deallocator>::type;
    static_assert(std::is_invocable_v<Deallocator, Container>,
                  "deallocator must have signature void(S::*)()");
    static_assert(
        std::is_same_v<std::invoke_result_t<Deallocator, Container>, void>,
        "deallocator must return void");

    struct CallbackImpl : public Callback {
      std::weak_ptr<Container> weak_ref;
      CallbackImpl(std::shared_ptr<Container> const &ref) : weak_ref{ref} {}
      ~CallbackImpl() override {
        if (auto const ref = weak_ref.lock()) {
          (ref.get()->*deallocator)();
        }
      }
    };

    static std::unique_ptr<Callback>
    make_callback(std::shared_ptr<Container> const &container) {
      return std::make_unique<CallbackImpl>(container);
    }
    friend class ResourceCleanup;
  };

  [[nodiscard]] auto size() const { return m_resources.size(); }
  [[nodiscard]] auto empty() const { return m_resources.empty(); }

  /**
   * @brief Register a resource for cleanup.
   * Internally, a weak pointer is created and stored in a callback.
   * @tparam Attorney   Specialization of @ref ResourceCleanup::Attorney
   *                    that wraps the class @c Container deallocator.
   * @tparam Container  Class that manages resources.
   */
  template <typename Attorney, typename Container>
  void push(std::shared_ptr<Container> const &resource) {
    m_resources.emplace_back(Attorney::make_callback(resource));
  }
};
