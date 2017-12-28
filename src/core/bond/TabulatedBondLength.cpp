#include "TabulatedBondLength.hpp"

/** Calculate a tabulated bond length force with number type_num (see
    \ref Bonded_ia_parameters) between particles p1 and p2 and add it
    to the particle forces. The force acts in the direction of the
    connecting vector between the particles. For distances smaller
    than the tabulated range it uses a linear extrapolation based on
    the first two tabulated force values.
    Needs feature TABULATED compiled in (see \ref config.hpp). */
int Bond::TabulatedBondLength::calc_bonded_pair_force(Particle *p1, Particle *p2, double dx[3], 
			  double force[3]) const 
{

  auto const dist = sqrt(sqrlen(dx));

  if (dist < m_tab_pot.cutoff()) {
    auto const fac = m_tab_pot.force(dist) / dist;

    for (int j = 0; j < 3; j++)
      force[j] -= fac * dx[j];

    return 0;
  } else {
    return 1;
  }

}

/** Calculate and return a tabulated bond length energy with number
    type_num (see \ref Bonded_ia_parameters) between particles p1 and
    p2. For distances smaller than the tabulated range it uses a
    quadratic extrapolation based on the first two tabulated force
    values and the first tabulated energy value.
    Needs feature TABULATED compiled in (see \ref config.hpp). */
int Bond::TabulatedBondLength::calc_bonded_pair_energy(Particle *p1, Particle *p2, 
			   double dx[3], double *_energy) const
{

  double dist = sqrt(sqrlen(dx));

  if (dist < m_tab_pot.cutoff()) {
    *_energy = m_tab_pot.energy(dist);
    return 0;
  } else {
    return 1;
  }

}
