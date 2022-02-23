#ifndef ESPRESSO_EKREACTION_HPP
#define ESPRESSO_EKREACTION_HPP

#include "EKReactant.hpp"
#include "EKinWalberlaBase.hpp"
#include "LatticeWalberla.hpp"

#include "generated_kernels/electrokinetics/reactions/ReactionKernel_1.h"
#include "generated_kernels/electrokinetics/reactions/ReactionKernel_2.h"
#include "generated_kernels/electrokinetics/reactions/ReactionKernel_3.h"
#include "generated_kernels/electrokinetics/reactions/ReactionKernel_4.h"
#include "generated_kernels/electrokinetics/reactions/ReactionKernel_5.h"

#include <memory>

namespace walberla {

namespace detail {
template <typename FloatType>
auto get_reaction_details(
    const std::shared_ptr<EKReactant<FloatType>> &reactant) {
  const auto order = reactant->get_order();
  const auto stoech_coeff = reactant->get_stoech_coeff();
  const auto density_id =
      walberla::BlockDataID(reactant->get_species()->get_density_id());

  return std::make_tuple(density_id, order, stoech_coeff);
}
} // namespace detail

template <typename FloatType> class EKReaction {
private:
  std::vector<std::shared_ptr<EKReactant<FloatType>>> m_reactants;
  FloatType m_coefficient;

  std::shared_ptr<LatticeWalberla> m_lattice;

public:
  EKReaction(std::shared_ptr<LatticeWalberla> lattice,
             std::vector<std::shared_ptr<EKReactant<FloatType>>> reactants,
             FloatType coefficient)
      : m_reactants(std::move(reactants)), m_coefficient(coefficient),
        m_lattice(std::move(lattice)) {}

  void set_coefficient(FloatType coefficient) { m_coefficient = coefficient; }
  [[nodiscard]] FloatType get_coefficient() const { return m_coefficient; }

  auto get_kernel() const {
    switch (m_reactants.size()) {
    case 1: {
      const auto &reactant = m_reactants[0];
      const auto [density_id_0, order_0, stoech_coeff_0] =
          detail::get_reaction_details<FloatType>(reactant);

      auto kernel = std::make_shared<pystencils::ReactionKernel_1>(
          density_id_0, order_0, get_coefficient(), stoech_coeff_0);

      return kernel->getSweep();
    }
    case 2: {
      const auto [density_id_0, order_0, stoech_coeff_0] =
          detail::get_reaction_details<FloatType>(m_reactants[0]);
      const auto [density_id_1, order_1, stoech_coeff_1] =
          detail::get_reaction_details<FloatType>(m_reactants[1]);

      auto kernel = std::make_shared<pystencils::ReactionKernel_2>(
          density_id_0, density_id_1, order_0, order_1, get_coefficient(),
          stoech_coeff_0, stoech_coeff_1);

      return kernel->getSweep();
    }
    case 3: {
      const auto [density_id_0, order_0, stoech_coeff_0] =
          detail::get_reaction_details<FloatType>(m_reactants[0]);
      const auto [density_id_1, order_1, stoech_coeff_1] =
          detail::get_reaction_details<FloatType>(m_reactants[1]);
      const auto [density_id_2, order_2, stoech_coeff_2] =
          detail::get_reaction_details<FloatType>(m_reactants[2]);

      auto kernel = std::make_shared<pystencils::ReactionKernel_3>(
          density_id_0, density_id_1, density_id_2, order_0, order_1, order_2,
          get_coefficient(), stoech_coeff_0, stoech_coeff_1, stoech_coeff_2);

      return kernel->getSweep();
    }
    case 4: {
      const auto [density_id_0, order_0, stoech_coeff_0] =
          detail::get_reaction_details<FloatType>(m_reactants[0]);
      const auto [density_id_1, order_1, stoech_coeff_1] =
          detail::get_reaction_details<FloatType>(m_reactants[1]);
      const auto [density_id_2, order_2, stoech_coeff_2] =
          detail::get_reaction_details<FloatType>(m_reactants[2]);
      const auto [density_id_3, order_3, stoech_coeff_3] =
          detail::get_reaction_details<FloatType>(m_reactants[3]);

      auto kernel = std::make_shared<pystencils::ReactionKernel_4>(
          density_id_0, density_id_1, density_id_2, density_id_3, order_0,
          order_1, order_2, order_3, get_coefficient(), stoech_coeff_0,
          stoech_coeff_1, stoech_coeff_2, stoech_coeff_3);

      return kernel->getSweep();
    }
    case 5: {
      const auto [density_id_0, order_0, stoech_coeff_0] =
          detail::get_reaction_details<FloatType>(m_reactants[0]);
      const auto [density_id_1, order_1, stoech_coeff_1] =
          detail::get_reaction_details<FloatType>(m_reactants[1]);
      const auto [density_id_2, order_2, stoech_coeff_2] =
          detail::get_reaction_details<FloatType>(m_reactants[2]);
      const auto [density_id_3, order_3, stoech_coeff_3] =
          detail::get_reaction_details<FloatType>(m_reactants[3]);
      const auto [density_id_4, order_4, stoech_coeff_4] =
          detail::get_reaction_details<FloatType>(m_reactants[4]);

      auto kernel = std::make_shared<pystencils::ReactionKernel_5>(
          density_id_0, density_id_1, density_id_2, density_id_3, density_id_4,
          order_0, order_1, order_2, order_3, order_4, get_coefficient(),
          stoech_coeff_0, stoech_coeff_1, stoech_coeff_2, stoech_coeff_3,
          stoech_coeff_4);

      return kernel->getSweep();
    }
    default:
      throw std::runtime_error("reactions of this size are not implemented!");
    }
  }

  void perform_reaction() const {
    // TODO: if my understanding is correct:
    //  the kernels need to either run in the ghost layers and do the
    //  synchronization before or not run and do a synchronization afterwards.
    //  The better solution is probably the latter one. Not sure why it fails
    //  atm.

    auto kernel = get_kernel();
    for (auto &block : *m_lattice->get_blocks()) {
      kernel(&block);
    }
  }
};
} // namespace walberla
#endif // ESPRESSO_EKREACTION_HPP
