#ifndef ESPRESSO_EK_REACTIONS_HPP
#define ESPRESSO_EK_REACTIONS_HPP

#include "EKReactionBase.hpp"
#include "EKReactions.hpp"

namespace EK {
extern EKReactions<walberla::EKReactionBase<double>> ek_reactions;

void perform_reactions();
} // namespace EK

#endif // ESPRESSO_EK_REACTIONS_HPP
