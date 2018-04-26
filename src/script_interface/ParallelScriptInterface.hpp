/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
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

#ifndef SCRIPT_INTERFACE_PARALLEL_SCRIPT_INTERFACE_HPP
#define SCRIPT_INTERFACE_PARALLEL_SCRIPT_INTERFACE_HPP

#include <utility>

#include "ParallelScriptInterfaceSlave.hpp"
#include "ScriptInterface.hpp"

namespace ScriptInterface {

class ParallelScriptInterface : public ScriptInterfaceBase {
  static Communication::MpiCallbacks *m_cb;

public:
  using CallbackAction = ParallelScriptInterfaceSlave::CallbackAction;

  ParallelScriptInterface(std::string const &name);
  ~ParallelScriptInterface() override;

  /**
   * @brief Initialize the mpi callback for instance creation.
   */
  static void initialize(Communication::MpiCallbacks &cb);

  bool operator==(ParallelScriptInterface const &rhs);
  bool operator!=(ParallelScriptInterface const &rhs);

  /**
   * @brief Get the payload object.
   */
  std::shared_ptr<ScriptInterfaceBase> get_underlying_object() const {
    return std::static_pointer_cast<ScriptInterfaceBase>(m_p);
  }

  void construct(VariantMap const &params) override;
  const std::string name() const { return m_p->name(); }
  void set_parameter(const std::string &name, const Variant &value) override;
  void
  set_parameters(const std::map<std::string, Variant> &parameters) override;
  ParameterMap valid_parameters() const override {
    return m_p->valid_parameters();
  }

  Variant get_parameter(const std::string &name) const override;
  VariantMap get_parameters() const override;
  Variant call_method(const std::string &name,
                      const VariantMap &parameters) override;

  /* Id mapping */
  Variant map_local_to_parallel_id(Variant const &value) const;
  Variant map_parallel_to_local_id(Variant const &value);

private:
  using map_t = std::map<ObjectId, std::shared_ptr<ParallelScriptInterface>>;

  VariantMap unwrap_variant_map(VariantMap const &map);

  void call(CallbackAction action, int has_params = 0) {
    m_cb->call(m_callback_id, static_cast<int>(action), has_params);
  }

  /**
   * @brief Remove instances that are not used by anybody but us.
   */
  void collect_garbage();

  /* Data members */
  int m_callback_id;
  /* Payload object */
  std::shared_ptr<ScriptInterfaceBase> m_p;
  map_t obj_map;
};

} /* namespace ScriptInterface */

#endif
