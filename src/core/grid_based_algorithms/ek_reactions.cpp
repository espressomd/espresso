#include "ek_reactions.hpp"

namespace EK {
EKReactions<walberla::EKReaction<double>> ek_reactions;

void perform_reactions() {
  if (ek_reactions.empty()) {
    return;
  }

  std::for_each(ek_reactions.begin(), ek_reactions.end(),
                [](auto const &reaction) { reaction->perform_reaction(); });
}
} // namespace EK
