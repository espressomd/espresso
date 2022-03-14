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

  [[nodiscard]] Utils::Vector3i
  get_mapped_index(const Utils::Vector3i &node) const {
    auto output = node;
    const auto shape = get_instance()->get_lattice()->get_grid_dimensions();
    for (auto i : {0, 1, 2}) {
      if (node[i] < 0) {
        output[i] = node[i] + shape[i];
      }

      if (output[i] < 0 or output[i] >= shape[i]) {
        auto constexpr formatter = Utils::Vector3i::formatter(", ");
        std::stringstream ss;
        ss << "provided index [" << formatter << node
           << "] is out of range for shape [" << formatter << shape << "]\n";
        throw std::runtime_error(ss.str());
      }
    }

    return output;
  }
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

  [[nodiscard]] Variant do_call_method(std::string const &method,
                                       VariantMap const &parameters) override {
    if (method == "add_node_to_boundary") {
      m_ekreaction->set_node_boundary(
          get_mapped_index(get_value<Utils::Vector3i>(parameters, "node")));
      return none;
    }
    if (method == "remove_node_from_boundary") {
      m_ekreaction->remove_node_from_boundary(
          get_mapped_index(get_value<Utils::Vector3i>(parameters, "node")));
      return none;
    }
    return none;
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
