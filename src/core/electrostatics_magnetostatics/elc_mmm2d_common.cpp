
#include "config.hpp"

#include "elc_mmm2d_common.hpp"
#include "utils.hpp"

#ifdef ELECTROSTATICS

double helper_first(int position, const double *partblk, const double *gblcblk, bool sum) {
  double sum1 = partblk[position + POQECM] * gblcblk[POQECP] +
                partblk[position + POQESM] * gblcblk[POQESP];

  double sum2 = partblk[position + POQECP] * gblcblk[POQECM] +
                partblk[position + POQESP] * gblcblk[POQESM];

  return sum ? sum1 + sum2 : sum1 - sum2;
}

double helper_second(int position, const double *partblk, const double *gblcblk, bool sum) {
  double sum1 = partblk[position + PQESSM] * gblcblk[PQESSP] +
                partblk[position + PQESCM] * gblcblk[PQESCP] +
                partblk[position + PQECSM] * gblcblk[PQECSP] +
                partblk[position + PQECCM] * gblcblk[PQECCP];

  double sum2 = partblk[position + PQESSP] * gblcblk[PQESSM] +
                partblk[position + PQESCP] * gblcblk[PQESCM] +
                partblk[position + PQECSP] * gblcblk[PQECSM] +
                partblk[position + PQECCP] * gblcblk[PQECCM];

  return sum ? sum1 + sum2 : sum1 - sum2;
}

#endif /* ELECTROSTATICS */
