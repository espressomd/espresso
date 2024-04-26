/*
 * Copyright (C) 2010-2022 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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

#ifndef COMMUNICATION_MPI_CALLBACKS
#define COMMUNICATION_MPI_CALLBACKS

/**
 * @file
 *
 * @ref Communication::MpiCallbacks manages MPI communication using a
 * visitor pattern. The program runs on the head node and is responsible
 * for calling callback functions on the worker nodes when necessary,
 * e.g. to broadcast global variables or run an algorithm in parallel.
 *
 * Callbacks are registered on the head node as function pointers via
 * the @ref REGISTER_CALLBACK. The visitor pattern allows using arbitrary
 * function signatures.
 */

#include <utils/NumeratedContainer.hpp>
#include <utils/tuple.hpp>
#include <utils/type_traits.hpp>

#include <boost/mpi/collectives/broadcast.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/packed_iarchive.hpp>
#include <boost/range/algorithm/remove_if.hpp>

#include <cassert>
#include <memory>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace Communication {

class MpiCallbacks;

namespace detail {
/**
 * @brief Check if a type can be used as a callback argument.
 *
 * This checks is a type can be a parameter type for a MPI callback.
 * Not allowed are pointers and non-const references, as output
 * parameters can not work across ranks.
 */
template <class T>
using is_allowed_argument =
    std::integral_constant<bool,
                           not(std::is_pointer_v<T> ||
                               (!std::is_const_v<std::remove_reference_t<T>> &&
                                std::is_lvalue_reference_v<T>))>;

template <class... Args>
using are_allowed_arguments =
    typename Utils::conjunction<is_allowed_argument<Args>...>::type;

/**
 * @brief Invoke a callable with arguments from an mpi buffer.
 *
 * @tparam F A Callable that can be called with Args as parameters.
 * @tparam Args Pack of arguments for @p F
 *
 * @param f Functor to be called
 * @param ia Buffer to extract the parameters from
 *
 * @return Return value of calling @p f.
 */
template <class F, class... Args>
auto invoke(F f, boost::mpi::packed_iarchive &ia) {
  static_assert(are_allowed_arguments<Args...>::value,
                "Pointers and non-const references are not allowed as "
                "arguments for callbacks.");

  /* This is the local receive buffer for the parameters. We have to strip
     away const so we can actually deserialize into it. */
  std::tuple<std::remove_const_t<std::remove_reference_t<Args>>...> params;
  Utils::for_each([&ia](auto &e) { ia >> e; }, params);

  /* We add const here, so that parameters can only be by value
     or const reference. Output parameters on callbacks are not
     sensible because the changes are not propagated back, so
     we make sure this does not compile. */
  return std::apply(f, std::as_const(params));
}

/**
 * @brief Type-erased interface for callbacks.
 *
 * This encapsulates the signature of the callback
 * and the parameter transfer, so that it can be
 * called without any type information on the parameters.
 */
struct callback_concept_t {
  /**
   * @brief Execute the callback.
   *
   * Unpack parameters for this callback, and then call it.
   */
  virtual void operator()(boost::mpi::communicator const &,
                          boost::mpi::packed_iarchive &) const = 0;
  virtual ~callback_concept_t() = default;
};

/**
 * @brief Callback without a return value.
 *
 * This is an implementation of a callback for a specific callable
 * @p F and a set of arguments to call it with.
 */
template <class F, class... Args>
struct callback_void_t final : public callback_concept_t {
  F m_f;

  callback_void_t(callback_void_t const &) = delete;
  callback_void_t(callback_void_t &&) = delete;

  template <class FRef>
  explicit callback_void_t(FRef &&f) : m_f(std::forward<FRef>(f)) {}
  void operator()(boost::mpi::communicator const &,
                  boost::mpi::packed_iarchive &ia) const override {
    detail::invoke<F, Args...>(m_f, ia);
  }
};

template <class F, class R, class... Args> struct FunctorTypes {
  using functor_type = F;
  using return_type = R;
  using argument_types = std::tuple<Args...>;
};

template <class C, class R, class... Args>
auto functor_types_impl(R (C::*)(Args...) const) {
  return FunctorTypes<C, R, Args...>{};
}

template <class F>
using functor_types =
    decltype(functor_types_impl(&std::remove_reference_t<F>::operator()));

template <class CRef, class C, class R, class... Args>
auto make_model_impl(CRef &&c, FunctorTypes<C, R, Args...>) {
  return std::make_unique<callback_void_t<C, Args...>>(std::forward<CRef>(c));
}

/**
 * @brief Make a @ref callback_model_t for a functor or lambda.
 *
 * The signature is deduced from F::operator() const, which has
 * to exist and can not be overloaded.
 */
template <typename F> auto make_model(F &&f) {
  return make_model_impl(std::forward<F>(f), functor_types<F>{});
}

/**
 * @brief Make a @ref callback_model_t for a function pointer.
 */
template <class... Args> auto make_model(void (*f_ptr)(Args...)) {
  return std::make_unique<callback_void_t<void (*)(Args...), Args...>>(f_ptr);
}
} // namespace detail

/**
 * @brief  The interface of the MPI callback mechanism.
 */
class MpiCallbacks {
public:
  /**
   * @brief RAII handle for a callback.
   *
   * This is what the client gets for registering a
   * dynamic (= not function pointer) callback.
   * It manages the lifetime of the callback handle
   * needed to call it. The handle has a type derived
   * from the signature of the callback, which makes
   * it possible to do static type checking on the
   * arguments.
   */
  template <class... Args> class CallbackHandle {
  public:
    template <typename F, class = std::enable_if_t<std::is_same_v<
                              typename detail::functor_types<F>::argument_types,
                              std::tuple<Args...>>>>
    CallbackHandle(std::shared_ptr<MpiCallbacks> cb, F &&f)
        : m_id(cb->add(std::forward<F>(f))), m_cb(std::move(cb)) {}

    CallbackHandle(CallbackHandle const &) = delete;
    CallbackHandle(CallbackHandle &&rhs) noexcept = default;
    CallbackHandle &operator=(CallbackHandle const &) = delete;
    CallbackHandle &operator=(CallbackHandle &&rhs) noexcept = default;

  private:
    int m_id;
    std::shared_ptr<MpiCallbacks> m_cb;

  public:
    /**
     * @brief Call the callback managed by this handle.
     *
     * The arguments are passed to the remote callees, it
     * must be possible to call the function with the provided
     * arguments, otherwise this will not compile.
     */
    template <class... ArgRef>
    auto operator()(ArgRef &&...args) const
        /* Enable if a hypothetical function with signature void(Args..)
         * could be called with the provided arguments. */
        -> std::enable_if_t<
            std::is_void_v<decltype(std::declval<void (*)(Args...)>()(
                std::forward<ArgRef>(args)...))>> {
      if (m_cb)
        m_cb->call(m_id, std::forward<ArgRef>(args)...);
    }

    ~CallbackHandle() {
      if (m_cb)
        m_cb->remove(m_id);
    }

    int id() const { return m_id; }
  };

  /* Avoid accidental copy, leads to mpi deadlock or split brain */
  MpiCallbacks(MpiCallbacks const &) = delete;
  MpiCallbacks &operator=(MpiCallbacks const &) = delete;

private:
  static auto &static_callbacks() {
    static std::vector<
        std::pair<void (*)(), std::unique_ptr<detail::callback_concept_t>>>
        callbacks;

    return callbacks;
  }

public:
  MpiCallbacks(boost::mpi::communicator comm,
               std::shared_ptr<boost::mpi::environment> mpi_env)
      : m_comm(std::move(comm)), m_mpi_env(std::move(mpi_env)) {
    /* Add a dummy at id 0 for loop abort. */
    m_callback_map.add(nullptr);

    for (auto &kv : static_callbacks()) {
      m_func_ptr_to_id[kv.first] = m_callback_map.add(kv.second.get());
    }
  }

  ~MpiCallbacks() {
    /* Release the clients on exit */
    if (m_comm.rank() == 0) {
      try {
        abort_loop();
      } catch (...) { // NOLINT(bugprone-empty-catch)
      }
    }
  }

private:
  /**
   * @brief Add a new callback.
   *
   * Add a new callback to the system. This is a collective
   * function that must be run on all nodes.
   *
   * @tparam F An object with a const call operator.
   *
   * @param f The callback function to add.
   * @return A handle with which the callback can be called.
   */
  template <typename F> auto add(F &&f) {
    m_callbacks.emplace_back(detail::make_model(std::forward<F>(f)));
    return m_callback_map.add(m_callbacks.back().get());
  }

public:
  /**
   * @brief Add a new callback.
   *
   * Add a new callback to the system. This is a collective
   * function that must be run on all nodes.
   *
   * @param fp Pointer to the static callback function to add.
   */
  template <class... Args> void add(void (*fp)(Args...)) {
    m_callbacks.emplace_back(detail::make_model(fp));
    const int id = m_callback_map.add(m_callbacks.back().get());
    m_func_ptr_to_id[reinterpret_cast<void (*)()>(fp)] = id;
  }

  /**
   * @brief Add a new callback.
   *
   * Add a new callback to the system. This is a collective
   * function that must be run on all nodes.
   *
   * @param fp Pointer to the static callback function to add.
   */
  template <class... Args> static void add_static(void (*fp)(Args...)) {
    static_callbacks().emplace_back(reinterpret_cast<void (*)()>(fp),
                                    detail::make_model(fp));
  }

private:
  /**
   * @brief Remove callback.
   *
   * Remove the callback id from the callback list.
   * This is a collective call that must be run on all nodes.
   *
   * @param id Identifier of the callback to remove.
   */
  void remove(int id) {
    m_callbacks.erase(
        boost::remove_if(m_callbacks,
                         [ptr = m_callback_map[id]](auto const &e) {
                           return e.get() == ptr;
                         }),
        m_callbacks.end());
    m_callback_map.remove(id);
  }

private:
  /**
   * @brief call a callback.
   *
   * Call the callback id.
   * The method can only be called on the head node
   * and has the prerequisite that the other nodes are
   * in the MPI loop.
   *
   * @param id The callback to call.
   * @param args Arguments for the callback.
   */
  template <class... Args> void call(int id, Args &&...args) const {
    if (m_comm.rank() != 0) {
      throw std::logic_error("Callbacks can only be invoked on rank 0.");
    }

    assert(m_callback_map.find(id) != m_callback_map.end() &&
           "m_callback_map and m_func_ptr_to_id disagree");

    /* Send request to worker nodes */
    boost::mpi::packed_oarchive oa(m_comm);
    oa << id;

    /* Pack the arguments into a packed mpi buffer. */
    Utils::for_each([&oa](auto &&e) { oa << e; },
                    std::forward_as_tuple(std::forward<Args>(args)...));

    boost::mpi::broadcast(m_comm, oa, 0);
  }

public:
  /**
   * @brief Call a callback on worker nodes.
   *
   * The callback is **not** called on the head node.
   *
   * This method can only be called on the head node.
   *
   * @param fp Pointer to the function to call.
   * @param args Arguments for the callback.
   */
  template <class... Args, class... ArgRef>
  auto call(void (*fp)(Args...), ArgRef &&...args) const ->
      /* enable only if fp can be called with the provided arguments */
      std::enable_if_t<std::is_void_v<decltype(fp(args...))>> {
    const int id = m_func_ptr_to_id.at(reinterpret_cast<void (*)()>(fp));

    call(id, std::forward<ArgRef>(args)...);
  }

  /**
   * @brief Call a callback on all nodes.
   *
   * This calls a callback on all nodes, including the head node.
   *
   * This method can only be called on the head node.
   *
   * @param fp Pointer to the function to call.
   * @param args Arguments for the callback.
   */
  template <class... Args, class... ArgRef>
  auto call_all(void (*fp)(Args...), ArgRef &&...args) const ->
      /* enable only if fp can be called with the provided arguments */
      std::enable_if_t<std::is_void_v<decltype(fp(args...))>> {
    call(fp, args...);
    fp(args...);
  }

  /**
   * @brief Start the MPI loop.
   *
   * This is the callback loop for the worker nodes. They block
   * on the MPI call and wait for a new callback request
   * coming from the head node.
   * This should be run on the worker nodes and must be running
   * so that the head node can issue call().
   */
  void loop() const {
    for (;;) {
      /* Communicate callback id and parameters */
      boost::mpi::packed_iarchive ia(m_comm);
      boost::mpi::broadcast(m_comm, ia, 0);

      int request;
      ia >> request;

      if (request == LOOP_ABORT) {
        break;
      }
      /* Call the callback */
      m_callback_map[request]->operator()(m_comm, ia);
    }
  }

  /**
   * @brief Abort the MPI loop.
   *
   * Make the worker nodes exit the MPI loop.
   */
  void abort_loop() { call(LOOP_ABORT); }

  /**
   * @brief The boost mpi communicator used by this instance
   */
  boost::mpi::communicator const &comm() const { return m_comm; }

  std::shared_ptr<boost::mpi::environment> share_mpi_env() const {
    return m_mpi_env;
  }

private:
  /**
   * @brief Id for the @ref abort_loop. Has to be 0.
   */
  static constexpr int LOOP_ABORT = 0;

  /**
   * The MPI communicator used for the callbacks.
   */
  boost::mpi::communicator m_comm;

  /**
   * The MPI environment used for the callbacks.
   */
  std::shared_ptr<boost::mpi::environment> m_mpi_env;

  /**
   * Internal storage for the callback functions.
   */
  std::vector<std::unique_ptr<detail::callback_concept_t>> m_callbacks;

  /**
   * Map of ids to callbacks.
   */
  Utils::NumeratedContainer<detail::callback_concept_t *> m_callback_map;

  /**
   * Mapping of function pointers to ids, so static callbacks can be
   * called by their pointer.
   */
  std::unordered_map<void (*)(), int> m_func_ptr_to_id;
};

template <class... Args>
using CallbackHandle = MpiCallbacks::CallbackHandle<Args...>;

/**
 * @brief Helper class to add callbacks before main.
 *
 * Should not be used directly, but via @ref REGISTER_CALLBACK.
 */
class RegisterCallback {

public:
  RegisterCallback() = delete;

  template <class... Args> explicit RegisterCallback(void (*cb)(Args...)) {
    MpiCallbacks::add_static(cb);
  }
};
} /* namespace Communication */

/**
 * @brief Register a static callback without return value.
 *
 * This registers a function as an mpi callback.
 * The macro should be used at global scope.
 *
 * @param cb A function
 */
#define REGISTER_CALLBACK(cb)                                                  \
  namespace Communication {                                                    \
  static ::Communication::RegisterCallback register_##cb(&(cb));               \
  }

#endif
