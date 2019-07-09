//
// Created by florian on 09.07.19.
//

#ifndef ESPRESSO_P3M_LOCALMESH_HPP
#define ESPRESSO_P3M_LOCALMESH_HPP

#include <utils/Vector.hpp>

/** Structure for local mesh parameters. */
struct LocalMesh {
  LocalMesh() = default;

  /**
   * @param my_left Lower corner of the local box
   * @param my_right Upper corner of the local box
   * @param halo Required halo size
   * @param ai Inverse grid constant
   * @param mesh_off Coordinates of the origin
   * @return
   */
  LocalMesh(const Utils::Vector3d &my_left, const Utils::Vector3d &my_right,
            const Utils::Vector3d &halo, const Utils::Vector3d &ai,
            const Utils::Vector3d &mesh_off) {
    /* inner left down grid point (global index) */
    for (int i = 0; i < 3; i++)
      in_ld[i] = (int)ceil(my_left[i] * ai[i] - mesh_off[i]);
    /* inner up right grid point (global index) */
    for (int i = 0; i < 3; i++)
      in_ur[i] = (int)floor(my_right[i] * ai[i] - mesh_off[i]);

    /* correct round-of errors at boundary */
    for (int i = 0; i < 3; i++) {
      if ((my_right[i] * ai[i] - mesh_off[i]) - in_ur[i] < ROUND_ERROR_PREC)
        in_ur[i]--;
      if (1.0 + (my_left[i] * ai[i] - mesh_off[i]) - in_ld[i] <
          ROUND_ERROR_PREC)
        in_ld[i]--;
    }
    /* inner grid dimensions */
    for (int i = 0; i < 3; i++)
      inner[i] = in_ur[i] - in_ld[i] + 1;
    /* index of left down grid point in global mesh */
    for (int i = 0; i < 3; i++)
      ld_ind[i] = (int)ceil((my_left[i] - halo[i]) * ai[i] - mesh_off[i]);
    /* left down margin */
    for (int i = 0; i < 3; i++)
      margin[i * 2] = in_ld[i] - ld_ind[i];
    /* up right grid point */
    int ind[3];
    for (int i = 0; i < 3; i++)
      ind[i] = (int)floor((my_right[i] + halo[i]) * ai[i] - mesh_off[i]);
    /* correct roundof errors at up right boundary */
    for (int i = 0; i < 3; i++)
      if (((my_right[i] + halo[i]) * ai[i] - mesh_off[i]) - ind[i] == 0)
        ind[i]--;
    /* up right margin */
    for (int i = 0; i < 3; i++)
      margin[(i * 2) + 1] = ind[i] - in_ur[i];

    /* grid dimension */
    size = 1;
    for (int i = 0; i < 3; i++) {
      dim[i] = ind[i] - ld_ind[i] + 1;
      size *= dim[i];
    }
    /* reduce inner grid indices from global to local */
    for (int i = 0; i < 3; i++)
      in_ld[i] = margin[i * 2];
    for (int i = 0; i < 3; i++)
      in_ur[i] = margin[i * 2] + inner[i];
  }

  /* local mesh characterization. */
  /** dimension (size) of local mesh. */
  Utils::Vector3i dim = {};
  /** number of local mesh points. */
  int size = 0;
  /** index of lower left corner of the
      local mesh in the global mesh. */
  Utils::Vector3i ld_ind = {};
  /** position of the first local mesh point. */
  Utils::Vector3d ld_pos = {};
  /** dimension of mesh inside node domain. */
  Utils::Vector3i inner = {};
  /** inner left down grid point */
  Utils::Vector3i in_ld = {};
  /** inner up right grid point + (1,1,1) */
  Utils::Vector3i in_ur = {};
  /** number of margin mesh points. */
  Utils::Array<int, 6> margin{};
};

#endif // ESPRESSO_LOCALMESH_HPP
