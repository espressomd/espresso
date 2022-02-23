#ifndef ESPRESSO_EK_REACTIONS_HPP
#define ESPRESSO_EK_REACTIONS_HPP

#include "EKReactions.hpp"
#include "electrokinetics/reactions/EKReactionBase.hpp"

namespace EK {
extern EKReactions<walberla::EKReactionBase<double>> ek_reactions;

void perform_reactions();
} // namespace EK

#endif // ESPRESSO_EK_REACTIONS_HPP
