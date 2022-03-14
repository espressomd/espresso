#ifndef SCRIPT_INTERFACE_WALBERLA_EKREACTION_HPP
#define SCRIPT_INTERFACE_WALBERLA_EKREACTION_HPP

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameter.hpp"

#include "EKReactant.hpp"
#include "LatticeWalberla.hpp"

#include "walberla_bridge/electrokinetics/reactions/EKReactionBase.hpp"
#include "walberla_bridge/electrokinetics/reactions/EKReactionImplBulk.hpp"
#include "walberla_bridge/electrokinetics/reactions/EKReactionIndexed.hpp"

#include <memory>

namespace ScriptInterface::walberla {
class EKReaction : public AutoParameters<::walberla::EKReactionBase<double>> {
public:
  [[nodiscard]] virtual std::shared_ptr<::walberla::EKReactionBase<double>>
  get_instance() const = 0;
};

class EKBulkReaction : public EKReaction {
public:
  void do_construct(VariantMap const &args) override {

    auto lattice =
        get_value<std::shared_ptr<LatticeWalberla>>(args, "lattice")->lattice();

    auto reactant = get_value<std::vector<Variant>>(args, "reactants");
    std::vector<std::shared_ptr<::walberla::EKReactant<double>>> output(
        reactant.size());
    auto get_instance = [](Variant &ekreactant) {
      return get_value<std::shared_ptr<EKReactant>>(ekreactant)->get_instance();
    };
    std::transform(reactant.begin(), reactant.end(), output.begin(),
                   get_instance);

    m_ekreaction = std::make_shared<::walberla::EKReactionImplBulk<double>>(
        lattice, output, get_value<double>(args, "coefficient"));

    add_parameters({{"coefficient",
                     [this](Variant const &v) {
                       m_ekreaction->set_coefficient(get_value<double>(v));
                     },
                     [this]() { return m_ekreaction->get_coefficient(); }}});
  }

  [[nodiscard]] std::shared_ptr<::walberla::EKReactionBase<double>>
  get_instance() const override {
    return m_ekreaction;
  }

private:
  std::shared_ptr<::walberla::EKReactionImplBulk<double>> m_ekreaction;
};

class EKIndexedReaction : public EKReaction {
public:
  void do_construct(VariantMap const &args) override {

    auto lattice =
        get_value<std::shared_ptr<LatticeWalberla>>(args, "lattice")->lattice();

    auto reactant = get_value<std::vector<Variant>>(args, "reactants");
    std::vector<std::shared_ptr<::walberla::EKReactant<double>>> output(
        reactant.size());
    auto get_instance = [](Variant &ekreactant) {
      return get_value<std::shared_ptr<EKReactant>>(ekreactant)->get_instance();
    };
    std::transform(reactant.begin(), reactant.end(), output.begin(),
                   get_instance);

    m_ekreaction = std::make_shared<::walberla::EKReactionIndexed<double>>(
        lattice, output, get_value<double>(args, "coefficient"));

    add_parameters({{"coefficient",
                     [this](Variant const &v) {
                       m_ekreaction->set_coefficient(get_value<double>(v));
                     },
                     [this]() { return m_ekreaction->get_coefficient(); }}});
  }

  [[nodiscard]] std::shared_ptr<::walberla::EKReactionBase<double>>
  get_instance() const override {
    return m_ekreaction;
  }

private:
  std::shared_ptr<::walberla::EKReactionIndexed<double>> m_ekreaction;
};
} // namespace ScriptInterface::walberla

#endif // SCRIPT_INTERFACE_WALBERLA_EKREACTION_HPP
