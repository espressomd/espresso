#include "Tabulated.hpp"
#include "utils.hpp" // for utils::malloc

/** Force factor lookup in a force table for bonded interactions (see
    \ref Bonded_ia_parameters). The force is calculated by linear
    interpolation between the closest tabulated values. There is no
    check for the upper bound! 
    Needs feature TABULATED compiled in (see \ref config.hpp).*/
double Bond::Tabulated::bonded_tab_force_lookup(double val) const
{
  int    ind;
  double dind;
  
  dind = (val - m_minval)*m_invstepsize;

  if( dind < 0.0 ) ind = 0;
  else ind = (int)dind;

  dind = dind - ind;
  /* linear interpolation between data points */
  return  m_f[ind]*(1.0-dind) + m_f[ind+1]*dind;
}


/** Energy lookup in a energy table for bonded interactions (see \ref
    Bonded_ia_parameters). The force is calculated by linear
    interpolation between the closest tabulated values. There is no
    check for the upper bound! 
    Needs feature TABULATED compiled in (see \ref config.hpp). */
double Bond::Tabulated::bonded_tab_energy_lookup(double val) const
{
  int ind;
  double dind;
  
  dind = (val - m_minval)*m_invstepsize;

  if( dind < 0.0 ) ind = 0;
  else ind = (int)dind;

  dind = dind - ind;
  /* linear interpolation between data points */
  return m_e[ind]*(1.0-dind) + m_e[ind+1]*dind;
}
