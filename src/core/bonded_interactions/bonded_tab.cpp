#include "bonded_interactions/bonded_tab.hpp"

#ifdef TABULATED
#include "communication.hpp"

int tabulated_bonded_set_params(int bond_type,
                                TabulatedBondedInteraction tab_type, double min,
                                double max, std::vector<double> const &energy,
                                std::vector<double> const &force) {
  if (bond_type < 0)
    return 1;

  assert(max >= min);
  assert((max == min) || force.size() > 1);
  assert(force.size() == energy.size());

  make_bond_type_exist(bond_type);

  /* set types */
  bonded_ia_params[bond_type].type = BONDED_IA_TABULATED;
  bonded_ia_params[bond_type].p.tab.type = tab_type;
  bonded_ia_params[bond_type].p.tab.pot = new TabulatedPotential;
  auto tab_pot = bonded_ia_params[bond_type].p.tab.pot;

  /* set number of interaction partners */
  if (tab_type == TAB_BOND_LENGTH) {
    tab_pot->minval = min;
    tab_pot->maxval = max;
    bonded_ia_params[bond_type].num = 1;
  } else if (tab_type == TAB_BOND_ANGLE) {
    tab_pot->minval = 0.0;
    tab_pot->maxval = PI + ROUND_ERROR_PREC;
    bonded_ia_params[bond_type].num = 2;
  } else if (tab_type == TAB_BOND_DIHEDRAL) {
    tab_pot->minval = 0.0;
    tab_pot->maxval = 2.0 * PI + ROUND_ERROR_PREC;
    bonded_ia_params[bond_type].num = 3;
  } else {
    runtimeError("Unsupported tabulated bond type.");
    return 1;
  }

  tab_pot->invstepsize = static_cast<double>(force.size() - 1) / (max - min);

  tab_pot->force_tab = force;
  tab_pot->energy_tab = energy;

  mpi_bcast_ia_params(bond_type, -1);

  return ES_OK;
}

#endif
