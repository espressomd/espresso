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
#include "utils/print.hpp"
#include "utils/tuple.hpp"

#include <functional>
#include <initializer_list>
#include <tuple>
#include <typeindex>

namespace Communication {
namespace detail {
template <class... Args> struct tuple {
  std::tuple<Args...> t;

  template <class Archive> void serialize(Archive &ar, const long int) {
    Utils::for_each([&ar](auto &e) { ar &e; }, t);
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
      Utils::print(__PRETTY_FUNCTION__);
    tuple<Args...> params;
    boost::mpi::broadcast(comm, params, 0);

    Utils::apply(m_f, params.t);
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
    template<class... Args>
    class CallbackHandle {
    public:
        CallbackHandle(CallbackHandle const&) = delete;
        CallbackHandle(CallbackHandle &&rhs) noexcept {
          std::swap(m_id, rhs.m_id);
          std::swap(m_cb, rhs.m_cb);
        }
        CallbackHandle& operator=(CallbackHandle const&) = delete;
        CallbackHandle& operator=(CallbackHandle &&rhs) noexcept {
          std::swap(m_id, rhs.m_id);
          std::swap(m_cb, rhs.m_cb);

          return *this;
        }

    private:
        int m_id;
        MpiCallbacks *m_cb;

        /* Enable if ArgRef pack matches Args up to reference. */
        template<class... ArgRef>
        using enable_if_args = std::enable_if_t<std::is_same<
                  std::tuple<std::remove_reference_t<ArgRef>...>,
                  std::tuple<Args...>
                >::value>;
    public:
        CallbackHandle() : m_id(0), m_cb(nullptr) {}
        CallbackHandle(int id, MpiCallbacks *cb) : m_id(id), m_cb(cb) {}

        template<class... ArgRef>
        enable_if_args<ArgRef...> operator()(ArgRef&&... args) const {
          if(m_cb) m_cb->call(m_id, std::forward<ArgRef>(args)...);
        }

        ~CallbackHandle() {
          if(m_cb) m_cb->remove(m_id);
        }
    };

    template<class F, class... Args>
    auto make_handle(int id, detail::model_t<F, Args...> const&) {
      return CallbackHandle<Args...>{id, this};
    }

  /* Avoid accidental copy, leads to mpi deadlock
     or split brain */
  MpiCallbacks(MpiCallbacks const &) = delete;
  MpiCallbacks &operator=(MpiCallbacks const &) = delete;

private:
    friend class RegisterCallback;
    static auto &static_callbacks() {
        static std::vector<std::pair<void(*)(), std::unique_ptr<detail::concept_t>>> m_callbacks;

        return m_callbacks;
    }

public:
  explicit MpiCallbacks(boost::mpi::communicator &comm,
                        bool abort_on_exit = true)
      : m_abort_on_exit(abort_on_exit), m_comm(comm) {
    /** Add a dummy at id 0 for loop abort. */
    m_callbacks.add(detail::make_model([](){}));

    for(auto& kv: static_callbacks()) {\
        m_func_ptr_to_id[kv.first] = m_callbacks.add(std::move(kv.second));
    }
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
  template <typename F> auto add(F &&f) {
    auto model = detail::make_model(std::forward<F>(f));
    return make_handle(m_callbacks.add(std::move(model)), *model);
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
  template <class... Args> void add(void (*fp)(Args...)) {
    const int id = m_callbacks.add(detail::make_model(fp));
    m_func_ptr_to_id[reinterpret_cast<void (*)()>(fp)] = id;
  }

  /**
   * @brief Remove callback.
   *
   * Remove the callback id from the callback list.
   * This is a collective call that must be run on all node.
   *
   * @param id Identifier of the callback to remove.
   */
  void remove(int id) {
      m_callbacks.remove(id);
  }

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
private:
  template <class... Args> void call(int id, Args &&... args) const {
      Utils::print(__PRETTY_FUNCTION__);
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

public:
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
  void loop() const {
      for (;;) {
          int request;
          /** Communicate callback id and parameters */
          boost::mpi::broadcast(m_comm, request, 0);
          /** id == 0 is loop_abort. */
          if (request == LOOP_ABORT) {
              break;
          } else {
              /** Call the callback */
              m_callbacks[request]->operator()(m_comm);
          }
      }
  }

  /**
   * @brief Abort mpi loop.
   *
   * Make the slaves exit the MPI loop.
   */
  void abort_loop() {
      call(LOOP_ABORT);
  }

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

    class RegisterCallback {

    public:
        RegisterCallback() = delete;

        template<class... Args>
        explicit RegisterCallback(void (*cb)(Args...)) {
            auto &cbs = MpiCallbacks::static_callbacks();

            cbs.emplace_back(reinterpret_cast<void(*)()>(cb), detail::make_model(cb));
        }
    };
} /* namespace Communication */


#define REGISTER_CALLBACK(cb) namespace Communication { static RegisterCallback register_##cb{&cb}; }

#endif
