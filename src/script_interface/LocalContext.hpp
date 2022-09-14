/*
 * Copyright (C) 2020-2022 The ESPResSo project
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
#ifndef SCRIPT_INTERFACE_LOCAL_CONTEXT_HPP
#define SCRIPT_INTERFACE_LOCAL_CONTEXT_HPP

#include "Context.hpp"
#include "ObjectHandle.hpp"
#include "ParallelExceptionHandler.hpp"

#include <utils/Factory.hpp>

#include <boost/mpi/communicator.hpp>

#include <cassert>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

namespace ScriptInterface {

/**
 * @brief Trivial context.
 *
 * This context just maintains a local copy of an object.
 */
class LocalContext : public Context {
  Utils::Factory<ObjectHandle> m_factory;
  bool m_is_head_node;
  boost::mpi::communicator const &m_comm;
  ParallelExceptionHandler m_parallel_exception_handler;

public:
  LocalContext(Utils::Factory<ObjectHandle> factory,
               boost::mpi::communicator const &comm)
      : m_factory(std::move(factory)), m_is_head_node(comm.rank() == 0),
        m_comm(comm),
        // NOLINTNEXTLINE(bugprone-throw-keyword-missing)
        m_parallel_exception_handler(comm) {}

  const Utils::Factory<ObjectHandle> &factory() const { return m_factory; }

  void notify_call_method(const ObjectHandle *, std::string const &,
                          VariantMap const &) override {}
  void notify_set_parameter(const ObjectHandle *, std::string const &,
                            Variant const &) override {}

  std::shared_ptr<ObjectHandle>
  make_shared(std::string const &name, const VariantMap &parameters) override {
    auto sp = m_factory.make(name);
    set_context(sp.get());

    sp->construct(parameters);

    return sp;
  }

  boost::string_ref name(const ObjectHandle *o) const override {
    assert(o);

    return factory().type_name(*o);
  }

  bool is_head_node() const override { return m_is_head_node; }
  void parallel_try_catch(std::function<void()> const &cb) const override {
    m_parallel_exception_handler.parallel_try_catch<std::exception>(cb);
  }
  boost::mpi::communicator const &get_comm() const override { return m_comm; }
};
} // namespace ScriptInterface

#endif
