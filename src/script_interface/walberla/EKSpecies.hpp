#ifndef SCRIPT_INTERFACE_WALBERLA__EKSPECIES_HPP
#define SCRIPT_INTERFACE_WALBERLA__EKSPECIES_HPP

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameter.hpp"

#include "grid_based_algorithms/walberla_blockforest.hpp"

#include "walberla_bridge/EKinWalberlaBase.hpp"
#include "walberla_bridge/ekin_walberla_init.hpp"

namespace ScriptInterface::walberla {

class EKSpecies : public AutoParameters<EKinWalberlaBase<double>> {
public:
  void do_construct(VariantMap const &args) override {
    m_ekinstance = ::walberla::new_ek_walberla(
        get_walberla_blockforest(), get_value<double>(args, "diffusion"),
        get_value<double>(args, "kT"), get_value<double>(args, "density"));

    add_parameters({{"diffusion",
                     [this](Variant const &v) {
                       m_ekinstance->set_diffusion(get_value<double>(v));
                     },
                     [this]() { return m_ekinstance->get_diffusion(); }},
                    {"kT",
                     [this](Variant const &v) {
                       m_ekinstance->set_kT(get_value<double>(v));
                     },
                     [this]() { return m_ekinstance->get_kT(); }}});
  }

  [[nodiscard]] std::shared_ptr<EKinWalberlaBase<double>> get_ekinstance() {
    return m_ekinstance;
  }

private:
  /* The actual constraint */
  std::shared_ptr<EKinWalberlaBase<double>> m_ekinstance;
};
} // namespace ScriptInterface::walberla

#endif // SCRIPT_INTERFACE_WALBERLA__EKSPECIES_HPP
