/*
  Copyright (C) 2010-2018 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
    Max-Planck-Institute for Polymer Research, Theory Group

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

#ifndef COMMUNICATION_MPI_CALLBACKS
#define COMMUNICATION_MPI_CALLBACKS

#include <boost/mpi/collectives/broadcast.hpp>
#include <boost/mpi/communicator.hpp>

#include "utils/NumeratedContainer.hpp"

#include <functional>
#include <initializer_list>
#include <tuple>

namespace Communication {

namespace detail {
namespace detail {
template <class F, class Tuple, std::size_t... I>
constexpr decltype(auto) apply_impl(F &&f, Tuple &&t,
                                    std::index_sequence<I...>) {
  return f(std::get<I>(std::forward<Tuple>(t))...);
}
} // namespace detail

template <class F, class Tuple>
constexpr decltype(auto) apply(F &&f, Tuple &&t) {
  return detail::apply_impl(
      std::forward<F>(f), std::forward<Tuple>(t),
      std::make_index_sequence<
          std::tuple_size<std::remove_reference_t<Tuple>>::value>{});
}

namespace detail {
template <class F, class Tuple, std::size_t... I>
constexpr void for_each_impl(F &&f, Tuple &t, std::index_sequence<I...>) {
  using expand = int[];
  (void)expand{0, ((void)(f(std::get<I>(t))), 0)...};
}
} // namespace detail

template <class F, class Tuple> void for_each(F &&f, Tuple &t) {
  detail::for_each_impl(
      std::forward<F>(f), t,
      std::make_index_sequence<std::tuple_size<Tuple>::value>{});
}

template <class... Args> struct tuple {
  std::tuple<Args...> t;

  template <class Archive> void serialize(Archive &ar, const long int) {
    for_each([&ar](auto &e) { ar &e; }, t);
  }
};

struct concept_t {
  virtual void operator()(const boost::mpi::communicator &) const = 0;

  virtual ~concept_t() = default;
};

template <class F, class... Args> struct model_t final : public concept_t {
  F m_f;

  template <class FRef>
  explicit model_t(FRef &&f) : m_f(std::forward<FRef>(f)) {}

  void operator()(const boost::mpi::communicator &comm) const override {
    tuple<Args...> params;
    boost::mpi::broadcast(comm, params, 0);

    apply(m_f, params.t);
  }
};

template <class CRef, class C, class R, class... Args>
auto make_model_impl(CRef &&c, R (C::*)(Args...) const) {
  return std::make_unique<model_t<C, Args...>>(std::forward<CRef>(c));
}

template <typename F> auto make_model(F &&f) {
  using C = std::remove_reference_t<F>;
  return make_model_impl(std::forward<F>(f), &C::operator());
}

template <class... Args> auto make_model(void (*f_ptr)(Args...)) {
  return std::make_unique<model_t<void (*)(Args...), Args...>>(f_ptr);
}
} // namespace detail

/**
 * @brief  The interface of the MPI callback mechanism.
 */
class MpiCallbacks {
public:
  /* Avoid accidental copy, leads to mpi deadlock
     or split brain */
  MpiCallbacks(MpiCallbacks const &) = delete;
  MpiCallbacks &operator=(MpiCallbacks const &) = delete;

public:
  explicit MpiCallbacks(boost::mpi::communicator &comm,
                        bool abort_on_exit = true)
      : m_abort_on_exit(abort_on_exit), m_comm(comm) {
    /** Add a dummy at id 0 for loop abort. */
    m_callbacks.add(detail::make_model([](){}));
  }

  ~MpiCallbacks() {
    /* Release the clients on exit */
    if (m_abort_on_exit && (m_comm.rank() == 0)) {
      abort_loop();
    }
  }

  /**
   * @brief Add a new callback.
   *
   * Add a new callback to the system. This is a collective
   * function that must be run on all nodes.
   *
   * @tparam F An object with a const call operator.
   *
   * @param f The callback function to add.
   * @return An integer id with which the callback can be called.
   **/
  template <typename F> int add(F &&f) {
    return m_callbacks.add(detail::make_model(std::forward<F>(f)));
  }

  /**
   * @brief Add a new static callback.
   *
   * Add a new callback to the system. This is a collective
   * function that must be run on all nodes.
   *
   * @param fp Pointer to the static callback function to add.
   * @return An integer id with which the callback can be called.
   **/
  template <class... Args> int add(void (*fp)(Args...)) {
    const int id = m_callbacks.add(detail::make_model(fp));
    m_func_ptr_to_id[reinterpret_cast<void (*)()>(fp)] = id;

    return id;
  }

  /**
   * @brief Remove callback.
   *
   * Remove the callback id from the callback list.
   * This is a collective call that must be run on all node.
   *
   * @param id Identifier of the callback to remove.
   */
  void remove(int id);

  /**
   * @brief call a callback.
   *
   * Call the callback id.
   * The method can only be called on the master
   * and has the prerequisite that the other nodes are
   * in the MPI loop.
   *
   * @param id The callback to call.
   * @param par1 First parameter to pass to the callback.
   * @param par2 Second parameter to pass to the callback.
   */
  template <class... Args> void call(int id, Args &&... args) const {
    /** Can only be call from master */
    if (m_comm.rank() != 0) {
      throw std::logic_error("Callbacks can only be invoked on rank 0.");
    }

    /** Check if callback exists */
    if (m_callbacks.find(id) == m_callbacks.end()) {
      throw std::out_of_range("Callback does not exists.");
    }

    /** Send request to slaves */
    boost::mpi::broadcast(m_comm, id, 0);
    detail::tuple<Args...> params{{args...}};
    boost::mpi::broadcast(m_comm, params, 0);
  }

  /**
   * @brief call a callback.
   *
   * Call a static callback by pointer.
   * The method can only be called the master
   * and has the prerequisite that the other nodes are
   * in the MPI loop.
   *
   * @param fp Static callback (e.g. the function name) to call.
   * @param par1 First parameter to pass to the callback.
   * @param par2 Second parameter to pass to the callback.
   */
  template <class... Args>
  void call(void (*fp)(std::decay_t<Args>...), Args &&... args) const {
    /** If the function pointer is invalid, map.at will throw
        an out_of_range exception. */
    const int id = m_func_ptr_to_id.at(reinterpret_cast<void (*)()>(fp));

    call(id, std::forward<Args>(args)...);
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
  void loop() const;

  /**
   * @brief Abort mpi loop.
   *
   * Make the slaves exit the MPI loop.
   */
  void abort_loop() const;

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
  Utils::NumeratedContainer<std::unique_ptr<detail::concept_t>> m_callbacks;

  /**
   * Mapping of function pointers to ids, so static callbacks can be
   * called by their pointer for backward compatibility.
   */
  std::unordered_map<void (*)(), int> m_func_ptr_to_id;
};
} /* namespace Communication */

#endif
