#include "EKReaction.hpp"

#include "EKinWalberlaBase.hpp"
#include "LatticeWalberla.hpp"

#include <memory>

namespace walberla {

template <typename FloatType>
std::function<void(IBlock *)> EKReaction<FloatType>::get_kernel() const {
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
  default:
    break;
  }

  throw std::runtime_error("reactions of this size are not implemented!");
}
} // namespace walberla
