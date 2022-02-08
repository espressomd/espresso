#ifndef ESPRESSO_EK_REACTIONS_HPP
#define ESPRESSO_EK_REACTIONS_HPP

#include "EKReaction.hpp"
#include "EKReactions.hpp"

namespace EK {
extern EKReactions<walberla::EKReaction<double>> ek_reactions;

void perform_reactions();
} // namespace EK

#endif // ESPRESSO_EK_REACTIONS_HPP
