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

#include <boost/mpi/communicator.hpp>
#include <boost/mpi/collectives/broadcast.hpp>

#include "utils/NumeratedContainer.hpp"
#include "utils/make_function.hpp"

#include <functional>
#include <initializer_list>
#include <tuple>

namespace Communication {

namespace detail {
    namespace detail {
        template <class F, class Tuple, std::size_t... I>
        constexpr decltype(auto) apply_impl(F&& f, Tuple&& t, std::index_sequence<I...>)
        {
            return f(std::get<I>(std::forward<Tuple>(t))...);
        }
    }  // namespace detail

    template <class F, class Tuple>
    constexpr decltype(auto) apply(F&& f, Tuple&& t)
    {
        return detail::apply_impl(
                std::forward<F>(f), std::forward<Tuple>(t),
                std::make_index_sequence<std::tuple_size<std::remove_reference_t<Tuple>>::value>{});
    }

    template<class... Args>
    struct tuple {
        std::tuple<Args...> t;

        template<class Archive, class Tuple, size_t... I>
        void serialize_impl(Archive &ar, Tuple &t, std::index_sequence<I...>) {
            int dummy[] = { 0, ((void)(ar & std::get<I>(t)), 0)... };
        }

        template<class Archive>
                void serialize(Archive &ar, const long int) {
                serialize_impl(ar, t, std::make_index_sequence<sizeof...(Args)>{});
                }
    };

    struct concept_t {
        virtual void operator()(const boost::mpi::communicator &) const = 0;

        virtual ~concept_t() = default;
    };

    template<class... Args>
    struct model_t final : public concept_t {
        using function_type = std::function<void(Args...)>;

        function_type m_f;

        explicit model_t(function_type f) : m_f(std::move(f)) {}

        void operator()(const boost::mpi::communicator &comm) const override {
          tuple<Args...> params;
          boost::mpi::broadcast(comm, params, 0);

          apply(m_f, params.t);
        }
    };

    template<class... Args>
    auto make_model(const std::function<void(Args...)> &f) {
        return std::make_unique<model_t<Args...>>(f);
    }
}

/**
 * @brief  The interface of the MPI callback mechanism.
 */
class MpiCallbacks {
  /* Avoid accidental copy, leads to mpi deadlock
     or split brain */
  MpiCallbacks(MpiCallbacks const &) = delete;
  MpiCallbacks &operator=(MpiCallbacks const &) = delete;

public:
  /** Function type of static callbacks. */
  typedef void (*func_ptr_type)(int, int);
  /** Type of the callback functions. */
  typedef std::function<void(int, int)> function_type;

  explicit MpiCallbacks(boost::mpi::communicator &comm,
                        bool abort_on_exit = true)
      : m_abort_on_exit(abort_on_exit), m_comm(comm) {
    /** Add a dummy at id 0 for loop abort. */
    m_callbacks.add(std::unique_ptr<detail::concept_t>{});
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
   * @param f The callback function to add.
   * @return An integer id with which the callback can be called.
   **/
  int add(const function_type &f);

  /**
   * @brief Add a new static callback.
   *
   * Add a new callback to the system. This is a collective
   * function that must be run on all nodes.
   *
   * @param fp Pointer to the static callback function to add.
   * @return An integer id with which the callback can be called.
   **/
  int add(func_ptr_type fp);

  /**
   * @brief Remove callback.
   *
   * Remove the callback id from the callback list.
   * This is a collective call that must be run on all node.
   *
   * @param id Identifier of the callback to remove.
   */
  void remove(const int id);

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
  void call(int id, int par1 = 0, int par2 = 0) const;

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
  void call(func_ptr_type fp, int par1 = 0, int par2 = 0) const;

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
  std::unordered_map<func_ptr_type, int> m_func_ptr_to_id;
};
} /* namespace Communication */

#endif
