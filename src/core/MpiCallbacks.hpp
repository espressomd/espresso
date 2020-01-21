/*
 * Copyright (C) 2010-2019 The ESPResSo project
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

#include <utils/NumeratedContainer.hpp>
#include <utils/as_const.hpp>
#include <utils/tuple.hpp>
#include <utils/type_traits.hpp>

#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/collectives/broadcast.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/range/algorithm/remove_if.hpp>

#include <functional>
#include <initializer_list>
#include <tuple>
#include <utility>

namespace Communication {
/**
 * This namespace contains tag types
 * to indicate what to do with the return
 * values of callbacks.
 */
namespace Result {
/** %Ignore result */
struct Ignore {};
constexpr auto ignore = Ignore{};
/** Return value from one rank */
struct OneRank {};
constexpr auto one_rank = OneRank{};
/** Return value from master rank */
struct MasterRank {};
constexpr auto master_rank = MasterRank{};
/** Reduce return value over all ranks */
struct Reduction {};
constexpr auto reduction = Reduction{};
} // namespace Result

namespace detail {
/**
 * @brief Check if a type can be used as a callback argument.
 *
 * This checks is a type can be a parameter type for a MPI callback.
 * Not allowed are pointers and non-const references, as output
 * parameters can not work across ranks.
 */
template <class T>
using is_allowed_argument = std::integral_constant<
    bool, not(std::is_pointer<T>::value ||
              (!std::is_const<std::remove_reference_t<T>>::value &&
               std::is_lvalue_reference<T>::value))>;

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
  return Utils::apply(f, Utils::as_const(params));
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
 * @p F and a set of arguments to call it with, where the there
 * is no return value.
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

/**
 * @brief Callback where the return value is ignored.
 *
 * This is an implementation of a callback for a specific callable
 * @p F and a set of arguments to call it with, where the valued from
 * all ranks are ignored.
 */
template <class F, class... Args>
struct callback_ignore_t final : public callback_concept_t {
  F m_f;

  callback_ignore_t(callback_ignore_t const &) = delete;
  callback_ignore_t(callback_ignore_t &&) = delete;

  template <class FRef>
  explicit callback_ignore_t(FRef &&f) : m_f(std::forward<FRef>(f)) {}
  void operator()(boost::mpi::communicator const &comm,
                  boost::mpi::packed_iarchive &ia) const override {
    detail::invoke<F, Args...>(m_f, ia);
  }
};

/**
 * @brief Callback with a return value from one rank.
 *
 * This is an implementation of a callback for a specific callable
 * @p F and a set of arguments to call it with, where the value from
 * one rank is returned. Only one node is allowed to return a value.
 */
template <class F, class... Args>
struct callback_one_rank_t final : public callback_concept_t {
  F m_f;

  callback_one_rank_t(callback_one_rank_t const &) = delete;
  callback_one_rank_t(callback_one_rank_t &&) = delete;

  template <class FRef>
  explicit callback_one_rank_t(FRef &&f) : m_f(std::forward<FRef>(f)) {}
  void operator()(boost::mpi::communicator const &comm,
                  boost::mpi::packed_iarchive &ia) const override {
    auto const result = detail::invoke<F, Args...>(m_f, ia);

    assert(1 == boost::mpi::all_reduce(comm, static_cast<int>(!!result),
                                       std::plus<>()) &&
           "Incorrect number of return values");

    /* If this rank returned a result, send it to the head node. */
    if (!!result) {
      comm.send(0, 42, *result);
    }
  }
};

/**
 * @brief Callback with a return value from the master rank.
 *
 * This is an implementation of a callback for a specific callable
 * @p F and a set of arguments to call it with, where the value from
 * the master rank is returned.
 */
template <class F, class... Args>
struct callback_master_rank_t final : public callback_concept_t {
  F m_f;

  callback_master_rank_t(callback_master_rank_t const &) = delete;
  callback_master_rank_t(callback_master_rank_t &&) = delete;

  template <class FRef>
  explicit callback_master_rank_t(FRef &&f) : m_f(std::forward<FRef>(f)) {}
  void operator()(boost::mpi::communicator const &comm,
                  boost::mpi::packed_iarchive &ia) const override {
    (void)detail::invoke<F, Args...>(m_f, ia);
  }
};

/**
 * @brief Callback with return value reduction.
 *
 * This is an implementation of a callback for a specific callable
 * @p F and a set of arguments to call it with, where the return
 * value is reduced over the communicator.
 */
template <class Op, class F, class... Args>
struct callback_reduce_t final : public callback_concept_t {
  Op m_op;
  F m_f;

  template <class OpRef, class FRef>
  explicit callback_reduce_t(OpRef &&op, FRef &&f)
      : m_op(std::forward<OpRef>(op)), m_f(std::forward<FRef>(f)) {}

  /**
   * @brief Execute the callback.
   *
   * Receive parameters for this callback, and then call it.
   *
   * @param comm The communicator to receive the parameters on.
   */
  void operator()(boost::mpi::communicator const &comm,
                  boost::mpi::packed_iarchive &ia) const override {
    /* Call the callback function, and reduce over the results with
     * the stored reduction operation. */
    boost::mpi::reduce(comm, detail::invoke<F, Args...>(m_f, ia), m_op, 0);
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
auto make_model_impl(Result::Ignore, CRef &&c, FunctorTypes<C, R, Args...>) {
  return std::make_unique<callback_void_t<C, Args...>>(std::forward<CRef>(c));
}

/**
 * @brief Make a @ref callback_model_t for a functor or lambda.
 *
 * The signature is deduced from F::operator() const, which has
 * to exist and can not be overloaded.
 */
template <typename F> auto make_model(F &&f) {
  return make_model_impl(Result::Ignore{}, std::forward<F>(f),
                         functor_types<F>{});
}

/**
 * @brief Make a @ref callback_model_t for a function pointer.
 *
 * This instantiates an implementation of a callback for a function
 * pointer. The main task here is to transfer the signature from
 * the pointer to the callback_model_t by template argument type
 * deduction.
 */
template <class... Args> auto make_model(void (*f_ptr)(Args...)) {
  return std::make_unique<callback_void_t<void (*)(Args...), Args...>>(f_ptr);
}

template <class Op, class R, class... Args>
auto make_model(Result::Reduction, R (*f_ptr)(Args...), Op &&op) {
  return std::make_unique<
      callback_reduce_t<std::remove_reference_t<Op>, R (*)(Args...), Args...>>(
      std::forward<Op>(op), f_ptr);
}

template <class R, class... Args>
auto make_model(Result::Ignore, R (*f_ptr)(Args...)) {
  return std::make_unique<callback_ignore_t<R (*)(Args...), Args...>>(f_ptr);
}

template <class R, class... Args>
auto make_model(Result::OneRank, R (*f_ptr)(Args...)) {
  return std::make_unique<callback_one_rank_t<R (*)(Args...), Args...>>(f_ptr);
}

template <class R, class... Args>
auto make_model(Result::MasterRank, R (*f_ptr)(Args...)) {
  return std::make_unique<callback_master_rank_t<R (*)(Args...), Args...>>(
      f_ptr);
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
   * arguments. */
  template <class... Args> class CallbackHandle {
  public:
    template <typename F, class = std::enable_if_t<std::is_same<
                              typename detail::functor_types<F>::argument_types,
                              std::tuple<Args...>>::value>>
    CallbackHandle(MpiCallbacks *cb, F &&f)
        : m_id(cb->add(std::forward<F>(f))), m_cb(cb) {}

    CallbackHandle(CallbackHandle const &) = delete;
    CallbackHandle(CallbackHandle &&rhs) noexcept = default;
    CallbackHandle &operator=(CallbackHandle const &) = delete;
    CallbackHandle &operator=(CallbackHandle &&rhs) noexcept = default;

  private:
    int m_id;
    MpiCallbacks *m_cb;

  public:
    /**
     * @brief Call the callback managed by this handle.
     *
     * The arguments are passed to the remote callees, it
     * must be possible to call the function with the provided
     * arguments, otherwise this will not compile.
     */
    template <class... ArgRef>
    auto operator()(ArgRef &&... args) const
        /* Enable if a hypothetical function with signature void(Args..)
         * could be called with the provided arguments. */
        -> std::enable_if_t<
            std::is_void<decltype(std::declval<void (*)(Args...)>()(
                std::forward<ArgRef>(args)...))>::value> {
      assert(m_cb), m_cb->call(m_id, std::forward<ArgRef>(args)...);
    }

    ~CallbackHandle() {
      if (m_cb)
        m_cb->remove(m_id);
    }

    MpiCallbacks *cb() const { return m_cb; }
  };

  /* Avoid accidental copy, leads to mpi deadlock
     or split brain */
  MpiCallbacks(MpiCallbacks const &) = delete;
  MpiCallbacks &operator=(MpiCallbacks const &) = delete;

private:
  static auto &static_callbacks() {
    static std::vector<
        std::pair<void (*)(), std::unique_ptr<detail::callback_concept_t>>>
        m_callbacks;

    return m_callbacks;
  }

public:
  explicit MpiCallbacks(boost::mpi::communicator &comm,
                        bool abort_on_exit = true)
      : m_abort_on_exit(abort_on_exit), m_comm(comm) {
    /** Add a dummy at id 0 for loop abort. */
    m_callback_map.add(nullptr);

    for (auto &kv : static_callbacks()) {
      m_func_ptr_to_id[kv.first] = m_callback_map.add(kv.second.get());
    }
  }

  ~MpiCallbacks() {
    /* Release the clients on exit */
    if (m_abort_on_exit && (m_comm.rank() == 0)) {
      try {
        abort_loop();
      } catch (...) {
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

  /**
   * @brief Add a new callback with a return value.
   *
   * Add a new callback to the system. This is a collective
   * function that must be run on all nodes.
   * Tag is one of the tag types from @ref Communication::Result,
   * which indicates what to do with the return values.
   *
   * @param tag Tag type indicating return operation
   * @param tag_args Argument for the return operation, if any.
   * @param fp Pointer to the static callback function to add.
   */
  template <class Tag, class R, class... Args, class... TagArgs>
  static void add_static(Tag tag, R (*fp)(Args...), TagArgs &&... tag_args) {
    static_callbacks().emplace_back(
        reinterpret_cast<void (*)()>(fp),
        detail::make_model(tag, fp, std::forward<TagArgs>(tag_args)...));
  }

private:
  /**
   * @brief Remove callback.
   *
   * Remove the callback id from the callback list.
   * This is a collective call that must be run on all node.
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
   * The method can only be called on the master
   * and has the prerequisite that the other nodes are
   * in the MPI loop.
   *
   * @param id The callback to call.
   * @param args Arguments for the callback.
   */
  template <class... Args> void call(int id, Args &&... args) const {
    /* Can only be call from master */
    if (m_comm.rank() != 0) {
      throw std::logic_error("Callbacks can only be invoked on rank 0.");
    }

    /* Check if callback exists */
    if (m_callback_map.find(id) == m_callback_map.end()) {
      throw std::out_of_range("Callback does not exists.");
    }

    /* Send request to slaves */
    boost::mpi::packed_oarchive oa(m_comm);
    oa << id;

    /* Pack the arguments into a packed mpi buffer. */
    Utils::for_each([&oa](auto &&e) { oa << e; },
                    std::forward_as_tuple(std::forward<Args>(args)...));

    boost::mpi::broadcast(m_comm, oa, 0);
  }

public:
  /**
   * @brief call a callback.
   *
   * Call a static callback by pointer.
   * The method can only be called on the master
   * and has the prerequisite that the other nodes are
   * in the MPI loop. Also the function has to be previously
   * registered e.g. with the @ref REGISTER_CALLBACK macro.
   * The callback is not called on the head node.
   *
   * @param fp Pointer to the function to call.
   * @param args Arguments for the callback.
   */
  template <class... Args, class... ArgRef>
  auto call(void (*fp)(Args...), ArgRef &&... args) const ->
      /* Enable only if fp can be called with the provided arguments,
       * e.g. if fp(args...) is well-formed. */
      std::enable_if_t<std::is_void<decltype(fp(args...))>::value> {
    const int id = m_func_ptr_to_id.at(reinterpret_cast<void (*)()>(fp));

    call(id, std::forward<ArgRef>(args)...);
  }

  /**
   * @brief call a callback.
   *
   * Call a static callback by pointer.
   * The method can only be called on the master
   * and has the prerequisite that the other nodes are
   * in the MPI loop. Also the function has to be previously
   * registered e.g. with the @ref REGISTER_CALLBACK macro.
   * The callback is called on the head node.
   *
   * @param fp Pointer to the function to call.
   * @param args Arguments for the callback.
   */
  template <class... Args, class... ArgRef>
  void call_all(void (*fp)(Args...), ArgRef &&... args) const {
    call(fp, args...);
    fp(args...);
  }

  /**
   * @brief Call a callback and reduce the result over all nodes.
   *
   * This calls a callback on all nodes, including the head node,
   * and does a mpi reduction with the registered operation.
   *
   * This method can only be called on the head node.
   */
  template <class Op, class R, class... Args>
  auto call(Result::Reduction, Op op, R (*fp)(Args...), Args... args) const
      -> std::remove_reference_t<decltype(op(std::declval<R>(),
                                             std::declval<R>()))> {
    const int id = m_func_ptr_to_id.at(reinterpret_cast<void (*)()>(fp));

    call(id, args...);

    std::remove_cv_t<std::remove_reference_t<decltype(
        op(std::declval<R>(), std::declval<R>()))>>
        result{};
    boost::mpi::reduce(m_comm, fp(args...), result, op, 0);

    return result;
  }

  /**
   * @brief Call a callback and ignore the result over all nodes.
   *
   * This calls a callback on all nodes, including the head node,
   * and ignore all return values.
   *
   * This method can only be called on the head node.
   */
  template <class R, class... Args, class... ArgRef>
  auto call(Result::Ignore, R (*fp)(Args...), ArgRef... args) const
      -> std::remove_reference_t<R> {

    const int id = m_func_ptr_to_id.at(reinterpret_cast<void (*)()>(fp));
    call(id, args...);

    fp(std::forward<Args>(args)...);

    return {};
  }

  /**
   * @brief Call a callback and reduce the result over all nodes.
   *
   * This calls a callback on all nodes, including the head node,
   * and returns the value from the node which returned an engaged
   * optional. Only one node is allowed to return a value.
   *
   * This method can only be called on the head node.
   */
  template <class R, class... Args, class... ArgRef>
  auto call(Result::OneRank, boost::optional<R> (*fp)(Args...),
            ArgRef... args) const -> std::remove_reference_t<R> {

    const int id = m_func_ptr_to_id.at(reinterpret_cast<void (*)()>(fp));
    call(id, args...);

    auto const local_result = fp(std::forward<Args>(args)...);

    assert(1 == boost::mpi::all_reduce(m_comm, static_cast<int>(!!local_result),
                                       std::plus<>()) &&
           "Incorrect number of return values");

    if (!!local_result) {
      return *local_result;
    }

    std::remove_cv_t<std::remove_reference_t<R>> result;
    m_comm.recv(boost::mpi::any_source, boost::mpi::any_tag, result);
    return result;
  }

  /**
   * @brief Call a callback and return the result of the head node.
   *
   * This calls a callback on all nodes, including the head node,
   * and returns the value from the head node.
   *
   * This method can only be called on the head node.
   */
  template <class R, class... Args, class... ArgRef>
  auto call(Result::MasterRank, R (*fp)(Args...), ArgRef... args) const
      -> std::remove_reference_t<R> {

    const int id = m_func_ptr_to_id.at(reinterpret_cast<void (*)()>(fp));
    call(id, args...);

    return fp(std::forward<Args>(args)...);
  }

  /**
   * @brief Mpi slave loop.
   *
   * This is the callback loop for the slaves. They block
   * on the MPI call and wait for a new callback request
   * coming from the master.
   * This should be run on the slaves and must be running
   * so that the master can issue call().
   */
  void loop() const {
    for (;;) {
      /* Communicate callback id and parameters */
      boost::mpi::packed_iarchive ia(m_comm);
      boost::mpi::broadcast(m_comm, ia, 0);

      int request;
      ia >> request;

      /* id == 0 is loop_abort. */
      if (request == LOOP_ABORT) {
        break;
      }
      /* Call the callback */
      m_callback_map[request]->operator()(m_comm, ia);
    }
  }

  /**
   * @brief Abort mpi loop.
   *
   * Make the slaves exit the MPI loop.
   */
  void abort_loop() { call(LOOP_ABORT); }

  /**
   * @brief The boost mpi communicator used by this instance
   */
  boost::mpi::communicator const &comm() const { return m_comm; }

private:
  /**
   * @brief Id for the loop_abort. Has to be 0.
   */
  enum { LOOP_ABORT = 0 };

  /**
   * @brief If loop_abort should be called on destruction
   *        on the head node.
   */
  bool m_abort_on_exit;

  /**
   * The MPI communicator used for the callbacks.
   */
  boost::mpi::communicator &m_comm;

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

  template <class Tag, class R, class... Args, class... TagArgs>
  explicit RegisterCallback(Tag tag, R (*cb)(Args...), TagArgs &&... tag_args) {
    MpiCallbacks::add_static(tag, cb, std::forward<TagArgs>(tag_args)...);
  }
};
} /* namespace Communication */

/**
 * @brief Register a static callback without return value
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

/**
 * @brief Register a static callback whose return value is reduced
 *
 * This registers a function as an mpi callback with
 * reduction of the return values from the nodes.
 * The macro should be used at global scope.
 *
 * @param cb A function
 * @param op Reduction operation
 */
#define REGISTER_CALLBACK_REDUCTION(cb, op)                                    \
  namespace Communication {                                                    \
  static ::Communication::RegisterCallback                                     \
      register_reduction_##cb(::Communication::Result::Reduction{}, &(cb),     \
                              (op));                                           \
  }

/**
 * @brief Register a static callback whose return value is to be ignored
 *
 * This registers a function as an mpi callback with
 * ignored return values.
 * The macro should be used at global scope.
 *
 * @param cb A function
 */
#define REGISTER_CALLBACK_IGNORE(cb)                                           \
  namespace Communication {                                                    \
  static ::Communication::RegisterCallback                                     \
      register_ignore_##cb(::Communication::Result::Ignore{}, &(cb));          \
  }

/**
 * @brief Register a static callback which returns a value on only one node
 *
 * This registers a function as an mpi callback with
 * reduction of the return values from one node
 * where the value of the optional is set.
 * The macro should be used at global scope.
 *
 * @param cb A function
 */
#define REGISTER_CALLBACK_ONE_RANK(cb)                                         \
  namespace Communication {                                                    \
  static ::Communication::RegisterCallback                                     \
      register_one_rank_##cb(::Communication::Result::OneRank{}, &(cb));       \
  }

/**
 * @brief Register a static callback whose return value is ignored except on
 * the head node
 *
 * This registers a function as an mpi callback with
 * reduction of the return values from the head node.
 * The macro should be used at global scope.
 *
 * @param cb A function
 */
#define REGISTER_CALLBACK_MASTER_RANK(cb)                                      \
  namespace Communication {                                                    \
  static ::Communication::RegisterCallback                                     \
      register_master_rank_##cb(::Communication::Result::MasterRank{}, &(cb)); \
  }

#endif
