#ifndef SCRIPT_INTERFACE_WALBERLA_EKREACTANT_HPP
#define SCRIPT_INTERFACE_WALBERLA_EKREACTANT_HPP

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameter.hpp"

#include "EKSpecies.hpp"

#include "walberla_bridge/electrokinetics/reactions/EKReactant.hpp"

#include <memory>

namespace ScriptInterface::walberla {
class EKReactant : public AutoParameters<::walberla::EKReactant> {
public:
  void do_construct(VariantMap const &args) override {
    m_ekreactant = std::make_shared<::walberla::EKReactant>(
        get_value<std::shared_ptr<EKSpecies>>(args, "ekspecies")
            ->get_ekinstance(),
        get_value<double>(args, "stoech_coeff"),
        get_value<double>(args, "order"));

    add_parameters({{"order",
                     [this](Variant const &v) {
                       m_ekreactant->set_order(get_value<double>(v));
                     },
                     [this]() { return m_ekreactant->get_order(); }},
                    {"stoech_coeff",
                     [this](Variant const &v) {
                       m_ekreactant->set_stoech_coefficient(
                           get_value<double>(v));
                     },
                     [this]() { return m_ekreactant->get_stoech_coeff(); }}});
  }

  [[nodiscard]] std::shared_ptr<::walberla::EKReactant> get_instance() {
    return m_ekreactant;
  }

private:
  /* The actual instance */
  std::shared_ptr<::walberla::EKReactant> m_ekreactant;
};
} // namespace ScriptInterface::walberla

#endif // SCRIPT_INTERFACE_WALBERLA_EKREACTANT_HPP
