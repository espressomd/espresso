/*
  Copyright (C) 2010-2018 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
    Max-Planck-Institute for Polymer Research, Theory Group

  This file is part of ESPResSo.

  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/** \file
 *
 *  Routines, row decomposition, data structures and communication for the
 * 3D-FFT.
 *
 */
#include "fft-common.hpp"

#if defined(P3M) || defined(DP3M)

#include "communication.hpp"
#include "debug.hpp"

#include "utils.hpp"
#include "utils/memory.hpp"

#include <cstring>
#include <fftw3.h>

void fft_common_pre_init(fft_data_struct *fft) {
  for (int i = 0; i < 4; i++) {
    fft->plan[i].send_block = nullptr;
    fft->plan[i].send_size = nullptr;
    fft->plan[i].recv_block = nullptr;
    fft->plan[i].recv_size = nullptr;
  }

  fft->init_tag = 0;
  fft->max_comm_size = 0;
  fft->max_mesh_size = 0;
  fft->send_buf = nullptr;
  fft->recv_buf = nullptr;
  fft->data_buf = nullptr;
}

void fft_pack_block(double const *const in, double *const out,
                    int const start[3], int const size[3], int const dim[3],
                    int element) {
  /* mid and slow changing indices */
  int m, s;
  /* linear index of in grid, linear index of out grid */
  int li_in, li_out = 0;
  /* copy size */
  int copy_size;
  /* offsets for indices in input grid */
  int m_in_offset, s_in_offset;
  /* offsets for indices in output grid */
  int m_out_offset;

  copy_size = element * size[2] * sizeof(double);
  m_in_offset = element * dim[2];
  s_in_offset = element * (dim[2] * (dim[1] - size[1]));
  m_out_offset = element * size[2];
  li_in = element * (start[2] + dim[2] * (start[1] + dim[1] * start[0]));

  for (s = 0; s < size[0]; s++) {
    for (m = 0; m < size[1]; m++) {
      memmove(&(out[li_out]), &(in[li_in]), copy_size);
      li_in += m_in_offset;
      li_out += m_out_offset;
    }
    li_in += s_in_offset;
  }
}

void fft_pack_block_permute1(double const *const in, double *const out,
                             int const start[3], int const size[3],
                             int const dim[3], int element) {
  /* slow,mid and fast changing indices for input  grid */
  int s, m, f, e;
  /* linear index of in grid, linear index of out grid */
  int li_in, li_out = 0;
  /* offsets for indices in input grid */
  int m_in_offset, s_in_offset;
  /* offset for mid changing indices of output grid */
  int m_out_offset;

  m_in_offset = element * (dim[2] - size[2]);
  s_in_offset = element * (dim[2] * (dim[1] - size[1]));
  m_out_offset = (element * size[0]) - element;
  li_in = element * (start[2] + dim[2] * (start[1] + dim[1] * start[0]));

  for (s = 0; s < size[0]; s++) { /* fast changing out */
    li_out = element * s;
    for (m = 0; m < size[1]; m++) {   /* slow changing out */
      for (f = 0; f < size[2]; f++) { /* mid  changing out */
        for (e = 0; e < element; e++)
          out[li_out++] = in[li_in++];
        li_out += m_out_offset;
      }
      li_in += m_in_offset;
    }
    li_in += s_in_offset;
  }
}

void fft_pack_block_permute2(double const *const in, double *const out,
                             int const start[3], int const size[3],
                             int const dim[3], int element) {
  /* slow,mid and fast changing indices for input  grid */
  int s, m, f, e;
  /* linear index of in grid, linear index of out grid */
  int li_in, li_out = 0;
  /* offsets for indices in input grid */
  int m_in_offset, s_in_offset;
  /* offset for slow changing index of output grid */
  int s_out_offset;
  /* start index for mid changing index of output grid */
  int m_out_start;

  m_in_offset = element * (dim[2] - size[2]);
  s_in_offset = element * (dim[2] * (dim[1] - size[1]));
  s_out_offset = (element * size[0] * size[1]) - element;
  li_in = element * (start[2] + dim[2] * (start[1] + dim[1] * start[0]));

  for (s = 0; s < size[0]; s++) { /* mid changing out */
    m_out_start = element * (s * size[1]);
    for (m = 0; m < size[1]; m++) { /* fast changing out */
      li_out = m_out_start + element * m;
      for (f = 0; f < size[2]; f++) { /* slow  changing out */
        for (e = 0; e < element; e++)
          out[li_out++] = in[li_in++];
        li_out += s_out_offset;
      }
      li_in += m_in_offset;
    }
    li_in += s_in_offset;
  }
}

void fft_unpack_block(double const *const in, double *const out,
                      int const start[3], int const size[3], int const dim[3],
                      int element) {
  /* mid and slow changing indices */
  int m, s;
  /* linear index of in grid, linear index of out grid */
  int li_in = 0, li_out;
  /* copy size */
  int copy_size;
  /* offset for indices in input grid */
  int m_in_offset;
  /* offsets for indices in output grid */
  int m_out_offset, s_out_offset;

  copy_size = element * (size[2] * sizeof(double));
  m_out_offset = element * dim[2];
  s_out_offset = element * (dim[2] * (dim[1] - size[1]));
  m_in_offset = element * size[2];
  li_out = element * (start[2] + dim[2] * (start[1] + dim[1] * start[0]));

  for (s = 0; s < size[0]; s++) {
    for (m = 0; m < size[1]; m++) {
      memmove(&(out[li_out]), &(in[li_in]), copy_size);
      li_in += m_in_offset;
      li_out += m_out_offset;
    }
    li_out += s_out_offset;
  }
}

/************************************************
 * private functions
 ************************************************/

boost::optional<std::vector<int>> fft_find_comm_groups(const Vector3i &grid1,
                                                       const Vector3i &grid2,
                                                       const int *node_list1,
                                                       int *node_list2,
                                                       int *pos, int *my_pos) {
  int i;
  /* communication group cell size on grid1 and grid2 */
  int s1[3], s2[3];
  /* The communication group cells build the same super grid on grid1 and grid2
   */
  int ds[3];
  /* communication group size */
  int g_size = 1;
  /* comm. group cell index */
  int gi[3];
  /* position of a node in a grid */
  int p1[3], p2[3];
  /* node identity */
  int n;
  /* this_node position in the communication group. */
  int c_pos = -1;
  /* flag for group identification */
  int my_group = 0;

  /* calculate dimension of comm. group cells for both grids */
  if ((grid1[0] * grid1[1] * grid1[2]) != (grid2[0] * grid2[1] * grid2[2]))
    return boost::none; /* unlike number of nodes */
  for (i = 0; i < 3; i++) {
    s1[i] = grid1[i] / grid2[i];
    if (s1[i] == 0)
      s1[i] = 1;
    else if (grid1[i] != grid2[i] * s1[i])
      return boost::none; /* grids do not match!!! */

    s2[i] = grid2[i] / grid1[i];
    if (s2[i] == 0)
      s2[i] = 1;
    else if (grid2[i] != grid1[i] * s2[i])
      return boost::none; /* grids do not match!!! */

    ds[i] = grid2[i] / s2[i];
    g_size *= s2[i];
  }

  std::vector<int> group(g_size);

  /* calc node_list2 */
  /* loop through all comm. group cells */
  for (gi[2] = 0; gi[2] < ds[2]; gi[2]++)
    for (gi[1] = 0; gi[1] < ds[1]; gi[1]++)
      for (gi[0] = 0; gi[0] < ds[0]; gi[0]++) {
        /* loop through all nodes in that comm. group cell */
        for (i = 0; i < g_size; i++) {
          p1[0] = (gi[0] * s1[0]) + (i % s1[0]);
          p1[1] = (gi[1] * s1[1]) + ((i / s1[0]) % s1[1]);
          p1[2] = (gi[2] * s1[2]) + (i / (s1[0] * s1[1]));

          p2[0] = (gi[0] * s2[0]) + (i % s2[0]);
          p2[1] = (gi[1] * s2[1]) + ((i / s2[0]) % s2[1]);
          p2[2] = (gi[2] * s2[2]) + (i / (s2[0] * s2[1]));

          n = node_list1[get_linear_index(p1[0], p1[1], p1[2], grid1)];
          node_list2[get_linear_index(p2[0], p2[1], p2[2], grid2)] = n;

          pos[3 * n + 0] = p2[0];
          pos[3 * n + 1] = p2[1];
          pos[3 * n + 2] = p2[2];
          if (my_group == 1)
            group[i] = n;
          if (n == this_node && my_group == 0) {
            my_group = 1;
            c_pos = i;
            my_pos[0] = p2[0];
            my_pos[1] = p2[1];
            my_pos[2] = p2[2];
            i = -1; /* restart the loop */
          }
        }
        my_group = 0;
      }

  /* permute comm. group according to the nodes position in the group */
  /* This is necessary to have matching node pairs during communication! */
  while (c_pos > 0) {
    n = group[g_size - 1];
    for (i = g_size - 1; i > 0; i--)
      group[i] = group[i - 1];
    group[0] = n;
    c_pos--;
  }
  return group;
}

int fft_calc_local_mesh(int const n_pos[3], int const n_grid[3],
                        int const mesh[3], double const mesh_off[3],
                        int loc_mesh[3], int start[3]) {
  int i, last[3], size = 1;

  for (i = 0; i < 3; i++) {
    start[i] =
        (int)ceil((mesh[i] / (double)n_grid[i]) * n_pos[i] - mesh_off[i]);
    last[i] = (int)floor((mesh[i] / (double)n_grid[i]) * (n_pos[i] + 1) -
                         mesh_off[i]);
    /* correct round off errors */
    if ((mesh[i] / (double)n_grid[i]) * (n_pos[i] + 1) - mesh_off[i] - last[i] <
        1.0e-15)
      last[i]--;
    if (1.0 + (mesh[i] / (double)n_grid[i]) * n_pos[i] - mesh_off[i] -
            start[i] <
        1.0e-15)
      start[i]--;
    loc_mesh[i] = last[i] - start[i] + 1;
    size *= loc_mesh[i];
  }
  return size;
}

int fft_calc_send_block(int const pos1[3], int const grid1[3],
                        int const pos2[3], int const grid2[3],
                        int const mesh[3], double const mesh_off[3],
                        int block[6]) {
  int i, size = 1;
  int mesh1[3], first1[3], last1[3];
  int mesh2[3], first2[3], last2[3];

  fft_calc_local_mesh(pos1, grid1, mesh, mesh_off, mesh1, first1);
  fft_calc_local_mesh(pos2, grid2, mesh, mesh_off, mesh2, first2);

  for (i = 0; i < 3; i++) {
    last1[i] = first1[i] + mesh1[i] - 1;
    last2[i] = first2[i] + mesh2[i] - 1;
    block[i] = std::max(first1[i], first2[i]) - first1[i];
    block[i + 3] = (std::min(last1[i], last2[i]) - first1[i]) - block[i] + 1;
    size *= block[i + 3];
  }
  return size;
}

#endif /* defined(P3M) || defined(DP3M) */
