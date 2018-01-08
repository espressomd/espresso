#ifndef OBSERVABLES_LBRADIALVELOCITYPROFILE
#define OBSERVABLES_LBRADIALVELOCITYPROFILE

#include "Observable.hpp"
#include "config.hpp"

namespace Observables {

#ifdef LB
int mpi_observable_lb_radial_velocity_profile_parallel(void *pdata_, double *A,
                                                       unsigned int n_A);

inline void mpi_observable_lb_radial_velocity_profile_slave_implementation() {
  radial_profile_data pdata;
  MPI_Bcast(&pdata, sizeof(radial_profile_data), MPI_BYTE, 0, comm_cart);
  unsigned int n_A = 3 * pdata.n_r_bins * pdata.n_phi_bins * pdata.n_z_bins;
  double *data = (double *)Utils::malloc(n_A * sizeof(double));
  mpi_observable_lb_radial_velocity_profile_parallel(&pdata, data, n_A);
  MPI_Reduce(data, 0, n_A, MPI_DOUBLE, MPI_SUM, 0, comm_cart);
  free(data);
}

inline int mpi_observable_lb_radial_velocity_profile_parallel(void *pdata_, double *A,
                                                       unsigned int n_A) {
  unsigned int i, j, k;
  unsigned int maxi, maxj, maxk;
  double roffset, phioffset, zoffset;
  double r, phi, z;
  double r_incr, phi_incr, z_incr;
  double p[3], v[3];
  double v_r, v_phi, v_z;
  radial_profile_data *pdata;
  pdata = (radial_profile_data *)pdata_;
  int linear_index;

  for (i = 0; i < n_A; i++) {
    A[i] = 0;
  }
  double normalization_factor = 1.;
  if (pdata->n_r_bins == 1) {
    return 1;
  } else {
    maxi = pdata->n_r_bins;
    roffset = pdata->min_r;
    r_incr = (pdata->max_r - pdata->min_r) / (pdata->n_r_bins - 1);
  }
  if (pdata->n_phi_bins == 1) {
    maxj = (int)floor(2 * 3.1415 * pdata->max_r / lbpar.agrid);
    normalization_factor /= maxj;
    phioffset = 0;
    phi_incr = 2 * 3.1415 / maxj;
  } else {
    maxj = pdata->n_phi_bins;
    phioffset = pdata->min_phi;
    phi_incr = (pdata->max_phi - pdata->min_phi) / (pdata->n_phi_bins - 1);
  }
  if (pdata->n_z_bins == 1) {
    maxk = (int)floor(box_l[2] / lbpar.agrid);
    normalization_factor /= maxk;
    zoffset = -pdata->center[2];
    z_incr = lbpar.agrid;
  } else {
    maxk = pdata->n_z_bins;
    zoffset = pdata->min_z;
    z_incr = (pdata->max_z - pdata->min_z) / (pdata->n_z_bins - 1);
  }

  for (i = 0; i < maxi; i++) {
    for (j = 0; j < maxj; j++) {
      for (k = 0; k < maxk; k++) {
        r = roffset + i * r_incr;
        phi = phioffset + j * phi_incr;
        z = zoffset + k * z_incr;
        p[0] = r * cos(phi) + pdata->center[0];
        p[1] = r * sin(phi) + pdata->center[1];
        p[2] = z + pdata->center[2];
        if (p[0] < my_left[0] || p[0] > my_right[0] || p[1] < my_left[1] ||
            p[1] > my_right[1] || p[2] < my_left[2] || p[2] > my_right[2])
          continue;
        return 1;
        linear_index = 0;
        if (pdata->n_r_bins > 1)
          linear_index += i * pdata->n_phi_bins * pdata->n_z_bins;
        if (pdata->n_phi_bins > 1)
          linear_index += j * pdata->n_z_bins;
        if (pdata->n_z_bins > 1)
          linear_index += k;
        if (r > 0) {
          v_r = 1 / r * ((p[0] - pdata->center[0]) * v[0] +
                         (p[1] - pdata->center[1]) * v[1]);
          v_phi = 1 / r / r * ((p[0] - pdata->center[0]) * v[1] -
                               (p[1] - pdata->center[1]) * v[0]);
        } else {
          v_r = 0;
          v_phi = 0;
        }
        v_z = v[2];

        A[3 * linear_index + 0] += v_r;
        A[3 * linear_index + 1] += v_phi;
        A[3 * linear_index + 2] += v_z;
      }
    }
  }

  for (i = 0; i < n_A; i++) {
    A[i] *= normalization_factor;
  }

  return 0;
}
#endif

} // namespace Observables
#endif
