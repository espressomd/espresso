#ifndef ESPRESSO_ELC_MMM2D_COMMON_HPP
#define ESPRESSO_ELC_MMM2D_COMMON_HPP

#include "grid.hpp"

#ifdef ELECTROSTATICS

/** \name Product decomposition data organization
    For the cell blocks
    it is assumed that the lower blocks part is in the lower half.
    This has to have positive sign, so that has to be first. */
/*@{*/
#define POQESP 0
#define POQECP 1
#define POQESM 2
#define POQECM 3

#define PQESSP 0
#define PQESCP 1
#define PQECSP 2
#define PQECCP 3
#define PQESSM 4
#define PQESCM 5
#define PQECSM 6
#define PQECCM 7

/** inverse box dimensions */
/*@{*/
static double ux, ux2, uy, uy2, uz;
/*@}*/

/** struct to hold cached sin/cos values */
/*@{*/
typedef struct {
    double s, c;
} SCCache;

/*@}*/

static void elc_mmm2d_common_init_invBoxl() {
  ux = 1 / box_l[0];
  ux2 = ux * ux;
  uy = 1 / box_l[1];
  uy2 = uy * uy;
  uz = 1 / box_l[2];
}

inline double
elc_mmm2d_common_add_force_dir(int position, const double *partblk, const double *gblcblk) {
  return partblk[position + POQESM] * gblcblk[POQECP] -
         partblk[position + POQECM] * gblcblk[POQESP] +
         partblk[position + POQESP] * gblcblk[POQECM] -
         partblk[position + POQECP] * gblcblk[POQESM];
}

inline double
elc_mmm2d_common_add_force_z(int position, const double *partblk, const double *gblcblk) {
  return partblk[position + POQECM] * gblcblk[POQECP] +
         partblk[position + POQESM] * gblcblk[POQESP] -
         partblk[position + POQECP] * gblcblk[POQECM] -
         partblk[position + POQESP] * gblcblk[POQESM];
}

inline double
elc_mmm2d_common_dir_energy(int position, const double *partblk, const double *gblcblk) {
  return partblk[position + POQECM] * gblcblk[POQECP] +
         partblk[position + POQESM] * gblcblk[POQESP] +
         partblk[position + POQECP] * gblcblk[POQECM] +
         partblk[position + POQESP] * gblcblk[POQESM];
}

inline void
elc_mmm2d_common_PQ_setup(int position, int xCacheOffset, int yCacheOffset, double factor, double *partblk,
                          Utils::Span<const SCCache> scxcache, Utils::Span<const SCCache> scycache) {
  partblk[position + PQESSM] =
          scxcache[xCacheOffset].s * scycache[yCacheOffset].s / factor;
  partblk[position + PQESCM] =
          scxcache[xCacheOffset].s * scycache[yCacheOffset].c / factor;
  partblk[position + PQECSM] =
          scxcache[xCacheOffset].c * scycache[yCacheOffset].s / factor;
  partblk[position + PQECCM] =
          scxcache[xCacheOffset].c * scycache[yCacheOffset].c / factor;

  partblk[position + PQESSP] =
          scxcache[xCacheOffset].s * scycache[yCacheOffset].s * factor;
  partblk[position + PQESCP] =
          scxcache[xCacheOffset].s * scycache[yCacheOffset].c * factor;
  partblk[position + PQECSP] =
          scxcache[xCacheOffset].c * scycache[yCacheOffset].s * factor;
  partblk[position + PQECCP] =
          scxcache[xCacheOffset].c * scycache[yCacheOffset].c * factor;
}

inline double
elc_mmm2d_common_add_PQ_force_x(int position, const double *partblk, const double *gblcblk) {
  return partblk[position + PQESCM] * gblcblk[PQECCP] +
         partblk[position + PQESSM] * gblcblk[PQECSP] -
         partblk[position + PQECCM] * gblcblk[PQESCP] -
         partblk[position + PQECSM] * gblcblk[PQESSP] +
         partblk[position + PQESCP] * gblcblk[PQECCM] +
         partblk[position + PQESSP] * gblcblk[PQECSM] -
         partblk[position + PQECCP] * gblcblk[PQESCM] -
         partblk[position + PQECSP] * gblcblk[PQESSM];
}

inline double
elc_mmm2d_common_add_PQ_force_y(int position, const double *partblk, const double *gblcblk) {
  return partblk[position + PQECSM] * gblcblk[PQECCP] +
         partblk[position + PQESSM] * gblcblk[PQESCP] -
         partblk[position + PQECCM] * gblcblk[PQECSP] -
         partblk[position + PQESCM] * gblcblk[PQESSP] +
         partblk[position + PQECSP] * gblcblk[PQECCM] +
         partblk[position + PQESSP] * gblcblk[PQESCM] -
         partblk[position + PQECCP] * gblcblk[PQECSM] -
         partblk[position + PQESCP] * gblcblk[PQESSM];
}

inline double
elc_mmm2d_common_add_PQ_force_z(int position, const double *partblk, const double *gblcblk) {
  return partblk[position + PQECCM] * gblcblk[PQECCP] +
         partblk[position + PQECSM] * gblcblk[PQECSP] +
         partblk[position + PQESCM] * gblcblk[PQESCP] +
         partblk[position + PQESSM] * gblcblk[PQESSP] -
         partblk[position + PQECCP] * gblcblk[PQECCM] -
         partblk[position + PQECSP] * gblcblk[PQECSM] -
         partblk[position + PQESCP] * gblcblk[PQESCM] -
         partblk[position + PQESSP] * gblcblk[PQESSM];
}

inline double
elc_mmm2d_common_PQ_energy(int position, const double *partblk, const double *gblcblk) {
  return partblk[position + PQECCM] * gblcblk[PQECCP] +
         partblk[position + PQECSM] * gblcblk[PQECSP] +
         partblk[position + PQESCM] * gblcblk[PQESCP] +
         partblk[position + PQESSM] * gblcblk[PQESSP] +
         partblk[position + PQECCP] * gblcblk[PQECCM] +
         partblk[position + PQECSP] * gblcblk[PQECSM] +
         partblk[position + PQESCP] * gblcblk[PQESCM] +
         partblk[position + PQESSP] * gblcblk[PQESSM];
}

/** partblk[pos + #] = sccache.# * / factor alternating div/mul*/
/*@{*/
inline void
elc_mmm2d_common_setup(int position, double factor, int cacheOffset, double *partblk,
                       Utils::Span<const SCCache> sccache) {
  partblk[position + POQESM] = sccache[cacheOffset].s / factor;
  partblk[position + POQESP] = sccache[cacheOffset].s * factor;
  partblk[position + POQECM] = sccache[cacheOffset].c / factor;
  partblk[position + POQECP] = sccache[cacheOffset].c * factor;
}
/*@}*/

/* vector operations */

/** pdc = 0 */
inline void elc_mmm2d_common_clear_vec(double *pdc, int size) {
  for (int i = 0; i < size; i++)
    pdc[i] = 0;
}

/** pdc_d = pdc_s */
inline void elc_mmm2d_common_copy_vec(double *pdc_d, double *pdc_s, int size) {
  for (int i = 0; i < size; i++)
    pdc_d[i] = pdc_s[i];
}

/** pdc_d = pdc_s1 + pdc_s2 */
inline void elc_mmm2d_common_add_vec(double *pdc_d, double *pdc_s1, double *pdc_s2, int size) {
  for (int i = 0; i < size; i++)
    pdc_d[i] = pdc_s1[i] + pdc_s2[i];
}

/** pdc_d = scale*pdc_s1 + pdc_s2 */
inline void elc_mmm2d_common_addscale_vec(double *pdc_d, double scale, double *pdc_s1,
                                          double *pdc_s2, int size) {
  for (int i = 0; i < size; i++)
    pdc_d[i] = scale * pdc_s1[i] + pdc_s2[i];
}

/** pdc_d = scale*pdc */
inline void elc_mmm2d_common_scale_vec(double scale, double *pdc, int size) {
  for (int i = 0; i < size; i++)
    pdc[i] *= scale;
}

/* block indexing */

inline double *elc_mmm2d_common_block(double *p, int index, int size) {
  return &p[index * size];
}


#endif /* ELECTROSTATICS */
#endif /* ESPRESSO_ELC_MMM2D_COMMON_HPP */
