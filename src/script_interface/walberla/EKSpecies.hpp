#ifndef SCRIPT_INTERFACE_WALBERLA__EKSPECIES_HPP
#define SCRIPT_INTERFACE_WALBERLA__EKSPECIES_HPP

#include "LatticeWalberla.hpp"

#include "core/communication.hpp"
#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameter.hpp"

#include "walberla_bridge/EKinWalberlaBase.hpp"
#include "walberla_bridge/EKinWalberlaImpl.hpp"

namespace ScriptInterface::walberla {

template <typename T> T mpi_return_one_rank(const boost::optional<T> &val) {
  // on the head-rank
  if (comm_cart.rank() == 0) {
    // if value is already on head rank return it
    if (!!val)
      return *val;

    // otherwise receive it from the other ranks
    T result;
    comm_cart.recv(boost::mpi::any_source, 42, result);
    return result;
  }

  // if the result is on the other ranks send it to the head-rank
  if (!!val) {
    comm_cart.send(0, 42, *val);
  }
  // this is not necessary, only to silence "-Wreturn-type"
  return T{};
}

class EKSpecies : public AutoParameters<EKinWalberlaBase<double>> {
public:
  void do_construct(VariantMap const &args) override {
    m_lattice = get_value<std::shared_ptr<LatticeWalberla>>(args, "lattice");
    m_ekinstance = std::make_shared<::walberla::EKinWalberlaImpl<13, double>>(
        m_lattice->lattice(), get_value<double>(args, "diffusion"),
        get_value<double>(args, "kT"), get_value<double>(args, "valency"),
        get_value<Utils::Vector3d>(args, "ext_efield"),
        get_value<double>(args, "density"), get_value<bool>(args, "advection"),
        get_value<bool>(args, "friction_coupling"));

    add_parameters(
        {{"diffusion",
          [this](Variant const &v) {
            m_ekinstance->set_diffusion(get_value<double>(v));
          },
          [this]() { return m_ekinstance->get_diffusion(); }},
         {"kT",
          [this](Variant const &v) {
            m_ekinstance->set_kT(get_value<double>(v));
          },
          [this]() { return m_ekinstance->get_kT(); }},
         {"valency",
          [this](Variant const &v) {
            m_ekinstance->set_valency(get_value<double>(v));
          },
          [this]() { return m_ekinstance->get_valency(); }},
         {"ext_efield",
          [this](Variant const &v) {
            m_ekinstance->set_ext_efield(get_value<Utils::Vector3d>(v));
          },
          [this]() { return m_ekinstance->get_ext_efield(); }},
         {"advection",
          [this](Variant const &v) {
            m_ekinstance->set_advection(get_value<bool>(v));
          },
          [this]() { return m_ekinstance->get_advection(); }},
         {"friction_coupling",
          [this](Variant const &v) {
            m_ekinstance->set_friction_coupling(get_value<bool>(v));
          },
          [this]() { return m_ekinstance->get_friction_coupling(); }},
         {"shape", AutoParameter::read_only,
          [this]() {
            return m_ekinstance->get_lattice().get_grid_dimensions();
          }},
         {"lattice", AutoParameter::read_only,
          [this]() { return m_lattice; }}});
  }

  [[nodiscard]] std::shared_ptr<EKinWalberlaBase<double>> get_ekinstance() {
    return m_ekinstance;
  }

  [[nodiscard]] Variant do_call_method(std::string const &method,
                                       VariantMap const &parameters) override {
    if (method == "get_density") {
      return mpi_return_one_rank(m_ekinstance->get_node_density(
          get_value<Utils::Vector3i>(parameters, "position")));
    }
    if (method == "set_density") {
      m_ekinstance->set_node_density(
          get_value<Utils::Vector3i>(parameters, "position"),
          get_value<double>(parameters, "value"));
      m_ekinstance->ghost_communication();
      return none;
    }
    if (method == "is_boundary") {
      return mpi_return_one_rank(m_ekinstance->get_node_is_boundary(
          get_value<Utils::Vector3i>(parameters, "position"), false));
    }
    if (method == "set_node_flux_boundary") {
      m_ekinstance->set_node_flux_boundary(
          get_value<Utils::Vector3i>(parameters, "position"),
          get_value<Utils::Vector3d>(parameters, "flux"));
      return none;
    }
    if (method == "remove_node_flux_boundary") {
      m_ekinstance->remove_node_from_flux_boundary(
          get_value<Utils::Vector3i>(parameters, "position"));
      return none;
    }
    if (method == "get_node_flux_at_boundary") {
      const auto index = get_value<Utils::Vector3i>(parameters, "position");
      const auto result = m_ekinstance->get_node_is_flux_boundary(index);
      bool is_boundary = (result) ? *result : false;
      is_boundary =
          boost::mpi::all_reduce(comm_cart, is_boundary, std::logical_or<>());
      if (is_boundary) {
        return mpi_return_one_rank(
            m_ekinstance->get_node_flux_at_boundary(index));
      }
      return Variant{None{}};
    }
    if (method == "update_flux_boundary_from_shape") {
      m_ekinstance->update_flux_boundary_from_shape(
          get_value<std::vector<int>>(parameters, "raster_view"),
          get_value<std::vector<double>>(parameters, "value_view"));
      return none;
    }
    if (method == "update_density_boundary_from_shape") {
      m_ekinstance->update_density_boundary_from_shape(
          get_value<std::vector<int>>(parameters, "raster_view"),
          get_value<std::vector<double>>(parameters, "value_view"));
      return none;
    }
    if (method == "update_flux_boundary_from_list") {
      m_ekinstance->update_flux_boundary_from_list(
          get_value<std::vector<int>>(parameters, "nodes_view"),
          get_value<std::vector<double>>(parameters, "value_view"));
      return none;
    }
    if (method == "update_density_boundary_from_list") {
      m_ekinstance->update_density_boundary_from_list(
          get_value<std::vector<int>>(parameters, "nodes_view"),
          get_value<std::vector<double>>(parameters, "value_view"));
      return none;
    }
    return none;
  }

private:
  /* The actual constraint */
  std::shared_ptr<EKinWalberlaBase<double>> m_ekinstance;

  std::shared_ptr<LatticeWalberla> m_lattice;
};
} // namespace ScriptInterface::walberla

#endif // SCRIPT_INTERFACE_WALBERLA__EKSPECIES_HPP
