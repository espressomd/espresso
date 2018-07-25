#ifndef CORE_LB_MODES_FROM_POP_HPP
#define CORE_LB_MODES_FROM_POP_HPP

#include "Modes.hpp"
#include "PopulationView.hpp"

namespace LB {
template <class T, class MemoryLayout>
ESPRESSO_GPU_EXPORT Modes<T>
modes_from_pop(const PopulationView<T, MemoryLayout> &pop) {
  Modes<T> modes;

  auto const n0 = pop(0);
  auto const n1p = pop(1) + pop(2);
  auto const n1m = pop(1) - pop(2);
  auto const n2p = pop(3) + pop(4);
  auto const n2m = pop(3) - pop(4);
  auto const n3p = pop(5) + pop(6);
  auto const n3m = pop(5) - pop(6);
  auto const n4p = pop(7) + pop(8);
  auto const n4m = pop(7) - pop(8);
  auto const n5p = pop(9) + pop(10);
  auto const n5m = pop(9) - pop(10);
  auto const n6p = pop(11) + pop(12);
  auto const n6m = pop(11) - pop(12);
  auto const n7p = pop(13) + pop(14);
  auto const n7m = pop(13) - pop(14);
  auto const n8p = pop(15) + pop(16);
  auto const n8m = pop(15) - pop(16);
  auto const n9p = pop(17) + pop(18);
  auto const n9m = pop(17) - pop(18);

  /* mass mode */
  modes.mass = n0 + n1p + n2p + n3p + n4p + n5p + n6p + n7p + n8p + n9p;

  /* momentum modes */
  modes.momentum[0] = n1m + n4m + n5m + n6m + n7m;
  modes.momentum[1] = n2m + n4m - n5m + n8m + n9m;
  modes.momentum[2] = n3m + n6m - n7m + n8m - n9m;

  /* stress modes */
  modes.stress[0] = -n0 + n4p + n5p + n6p + n7p + n8p + n9p;
  modes.stress[1] = n1p - n2p + n6p + n7p - n8p - n9p;
  modes.stress[2] = n1p + n2p - n6p - n7p - n8p - n9p - 2. * (n3p - n4p - n5p);
  modes.stress[3] = n4p - n5p;
  modes.stress[4] = n6p - n7p;
  modes.stress[5] = n8p - n9p;

  /* kinetic modes */
  modes.kinetic[0] = -2. * n1m + n4m + n5m + n6m + n7m;
  modes.kinetic[1] = -2. * n2m + n4m - n5m + n8m + n9m;
  modes.kinetic[2] = -2. * n3m + n6m - n7m + n8m - n9m;
  modes.kinetic[3] = n4m + n5m - n6m - n7m;
  modes.kinetic[4] = n4m - n5m - n8m - n9m;
  modes.kinetic[5] = n6m - n7m - n8m + n9m;
  modes.kinetic[6] =
      n0 + n4p + n5p + n6p + n7p + n8p + n9p - 2. * (n1p + n2p + n3p);
  modes.kinetic[7] = -n1p + n2p + n6p + n7p - n8p - n9p;
  modes.kinetic[8] =
      -n1p - n2p - n6p - n7p - n8p - n9p + 2. * (n3p + n4p + n5p);

  return modes;
}
}

#endif
