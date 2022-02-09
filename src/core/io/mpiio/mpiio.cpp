/*
 * Copyright (C) 2010-2019 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
/** @file
 *
 * Concerning the file layouts.
 * - Scalar arrays are written like this:
 *   rank0 --- rank1 --- rank2 ...
 *   where each rank dumps its scalars in the ordering of the particles.
 * - Vector arrays are written in the rank ordering like scalar arrays.
 *   The ordering of the vector data is: v[0] v[1] v[2], so the data
 *   looks like this:
 *   v1[0] v1[1] v1[2] v2[0] v2[1] v2[2] v3[0] ...
 *
 * To be able to determine the rank boundaries (a multiple of
 * nlocalparts), the file 1.pref is written, which dumps the Exscan
 * results of nlocalparts, i.e. the prefixes in scalar arrays:
 * - 1.prefs looks like this:
 *   0 nlocalpats_rank0 nlocalparts_rank0+nlocalparts_rank1 ...
 *
 * Bonds are dumped as two arrays, namely 1.bond which stores the
 * bonding partners of the particles and 1.boff which stores the
 * iteration indices for each particle.
 * - 1.boff is a scalar array of size (nlocalpart + 1) per rank.
 * - The last element (at index nlocalpart) of 1.boff's subpart
 *   [rank * (nlocalpart + 1) : (rank + 1) * (nlocalpart + 1)]
 *   determines the number of bonds for processor "rank".
 * - In this subarray one can find the bonding partners of particle
 *   id[i]. The iteration indices for local part of 1.bonds are:
 *   subarray[i] : subarray[i+1]
 * - Take a look at the bond input code. It's easy to understand.
 */

#include "mpiio.hpp"

#include "Particle.hpp"
#include "bonded_interactions/bonded_interaction_data.hpp"
#include "cells.hpp"
#include "errorhandling.hpp"

#include <utils/Vector.hpp>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/stream.hpp>

#include <mpi.h>

#include <cerrno>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <string>
#include <sys/stat.h>
#include <tuple>
#include <utility>
#include <vector>

namespace Mpiio {

/**
 * @brief Fatal error handler.
 * On 1 MPI rank the error is recoverable and an exception is thrown.
 * On more than 1 MPI rank the error is not recoverable.
 */
static void fatal_error() {
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size == 1) {
    throw std::runtime_error("");
  }
  errexit();
}

/**
 * @brief Dump data @p arr of size @p len starting from prefix @p pref
 * of type @p T using @p MPI_T as MPI datatype. Beware, that @p T and
 * @p MPI_T have to match!
 *
 * @param fn The file name to write to (must not already exist!)
 * @param arr The array to dump
 * @param len The number of elements to dump
 * @param pref The prefix for this process
 * @param MPI_T The MPI datatype corresponding to the template parameter @p T
 */
template <typename T>
static void mpiio_dump_array(const std::string &fn, T *arr, std::size_t len,
                             std::size_t pref, MPI_Datatype MPI_T) {
  MPI_File f;
  int ret;

  ret = MPI_File_open(MPI_COMM_WORLD, const_cast<char *>(fn.c_str()),
                      // MPI_MODE_EXCL: Prohibit overwriting
                      MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_EXCL,
                      MPI_INFO_NULL, &f);
  if (ret) {
    char buf[MPI_MAX_ERROR_STRING];
    int buf_len;
    MPI_Error_string(ret, buf, &buf_len);
    buf[buf_len] = '\0';
    fprintf(stderr, "MPI-IO Error: Could not open file \"%s\": %s\n",
            fn.c_str(), buf);
    fatal_error();
  }
  ret = MPI_File_set_view(f, pref * sizeof(T), MPI_T, MPI_T,
                          const_cast<char *>("native"), MPI_INFO_NULL);
  ret |= MPI_File_write_all(f, arr, len, MPI_T, MPI_STATUS_IGNORE);
  MPI_File_close(&f);
  if (ret) {
    fprintf(stderr, "MPI-IO Error: Could not write file \"%s\".\n", fn.c_str());
    fatal_error();
  }
}

/**
 * @brief Dump the fields and bond information.
 * To be called by the head node only.
 *
 * @param fn The filename to write to
 * @param fields The dumped fields
 */
static void dump_info(const std::string &fn, unsigned fields) {
  FILE *f = fopen(fn.c_str(), "wb");
  if (!f) {
    fprintf(stderr, "MPI-IO Error: Could not open %s for writing.\n",
            fn.c_str());
    fatal_error();
  }
  static std::vector<int> npartners;
  int success = (fwrite(&fields, sizeof(fields), 1, f) == 1);
  // Pack the necessary information of bonded_ia_params:
  // The number of partners. This is needed to interpret the bond IntList.
  if (bonded_ia_params.size() > npartners.size())
    npartners.resize(bonded_ia_params.size());

  for (int i = 0; i < bonded_ia_params.size(); ++i) {
    npartners[i] = number_of_partners(*bonded_ia_params.at(i));
  }
  auto ia_params_size = static_cast<std::size_t>(bonded_ia_params.size());
  success =
      success && (fwrite(&ia_params_size, sizeof(std::size_t), 1, f) == 1);
  success =
      success && (fwrite(npartners.data(), sizeof(int), bonded_ia_params.size(),
                         f) == bonded_ia_params.size());
  fclose(f);
  if (!success) {
    fprintf(stderr, "MPI-IO Error: Failed to write %s.\n", fn.c_str());
    fatal_error();
  }
}

void mpi_mpiio_common_write(const std::string &prefix, unsigned fields,
                            const ParticleRange &particles) {
  auto const nlocalpart = static_cast<int>(particles.size());
  // Keep static buffers in order not having to allocate them on every
  // function call
  static std::vector<double> pos, vel;
  static std::vector<int> id, type;

  // Nlocalpart prefixes
  // Prefixes based for arrays: 3 * pref for vel, pos.
  int pref = 0, bpref = 0;
  MPI_Exscan(&nlocalpart, &pref, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  // Realloc static buffers if necessary
  if (nlocalpart > id.size())
    id.resize(nlocalpart);
  if (fields & MPIIO_OUT_POS && 3 * nlocalpart > pos.size())
    pos.resize(3 * nlocalpart);
  if (fields & MPIIO_OUT_VEL && 3 * nlocalpart > vel.size())
    vel.resize(3 * nlocalpart);
  if (fields & MPIIO_OUT_TYP && nlocalpart > type.size())
    type.resize(nlocalpart);

  // Pack the necessary information
  // Esp. rescale the velocities.
  int i1 = 0, i3 = 0;
  for (auto const &p : particles) {
    id[i1] = p.p.identity;
    if (fields & MPIIO_OUT_POS) {
      pos[i3] = p.r.p[0];
      pos[i3 + 1] = p.r.p[1];
      pos[i3 + 2] = p.r.p[2];
    }
    if (fields & MPIIO_OUT_VEL) {
      vel[i3] = p.m.v[0];
      vel[i3 + 1] = p.m.v[1];
      vel[i3 + 2] = p.m.v[2];
    }
    if (fields & MPIIO_OUT_TYP) {
      type[i1] = p.p.type;
    }
    i1++;
    i3 += 3;
  }

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0)
    dump_info(prefix + ".head", fields);
  mpiio_dump_array<int>(prefix + ".pref", &pref, 1, rank, MPI_INT);
  mpiio_dump_array<int>(prefix + ".id", id.data(), nlocalpart, pref, MPI_INT);
  if (fields & MPIIO_OUT_POS)
    mpiio_dump_array<double>(prefix + ".pos", pos.data(), 3 * nlocalpart,
                             3 * pref, MPI_DOUBLE);
  if (fields & MPIIO_OUT_VEL)
    mpiio_dump_array<double>(prefix + ".vel", vel.data(), 3 * nlocalpart,
                             3 * pref, MPI_DOUBLE);
  if (fields & MPIIO_OUT_TYP)
    mpiio_dump_array<int>(prefix + ".type", type.data(), nlocalpart, pref,
                          MPI_INT);

  if (fields & MPIIO_OUT_BND) {
    std::vector<char> bonds;

    /* Construct archive that pushes back to the bond buffer */
    {
      namespace io = boost::iostreams;
      io::stream_buffer<io::back_insert_device<std::vector<char>>> os{
          io::back_inserter(bonds)};
      boost::archive::binary_oarchive bond_archiver{os};

      for (auto const &p : particles) {
        bond_archiver << p.bonds();
      }
    }

    // Determine the prefixes in the bond file
    int bonds_size = static_cast<int>(bonds.size());
    MPI_Exscan(&bonds_size, &bpref, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    mpiio_dump_array<int>(prefix + ".boff", &bonds_size, 1, rank, MPI_INT);
    mpiio_dump_array<char>(prefix + ".bond", bonds.data(), bonds.size(), bpref,
                           MPI_CHAR);
  }
}

/**
 * @brief Get the number of elements in a file by its file size and @p elem_sz.
 * I.e. query the file size using stat(2) and divide it by @p elem_sz.
 *
 * @param fn The filename
 * @param elem_sz Size of a single element
 * @return The number of elements stored in the file
 */
static int get_num_elem(const std::string &fn, std::size_t elem_sz) {
  // Could also be done via MPI_File_open, MPI_File_get_size,
  // MPI_File_close.
  struct stat st;
  errno = 0;
  if (stat(fn.c_str(), &st) != 0) {
    fprintf(stderr, "MPI-IO Input Error: Could not get file size of %s: %s\n",
            fn.c_str(), strerror(errno));
    fatal_error();
  }
  return static_cast<int>(st.st_size / elem_sz);
}

/**
 * @brief Read a previously dumped array of size @p len starting from prefix
 * @p pref of type @p T using @p MPI_T as MPI datatype. Beware, that
 * @p T and @p MPI_T have to match!
 *
 * @param fn The file name to read from
 * @param arr The array to populate
 * @param len The number of elements to read
 * @param pref The prefix for this process
 * @param MPI_T The MPI datatype corresponding to the template parameter @p T
 */
template <typename T>
static void mpiio_read_array(const std::string &fn, T *arr, std::size_t len,
                             std::size_t pref, MPI_Datatype MPI_T) {
  MPI_File f;
  int ret;

  ret = MPI_File_open(MPI_COMM_WORLD, const_cast<char *>(fn.c_str()),
                      MPI_MODE_RDONLY, MPI_INFO_NULL, &f);

  if (ret) {
    char buf[MPI_MAX_ERROR_STRING];
    int buf_len;
    MPI_Error_string(ret, buf, &buf_len);
    buf[buf_len] = '\0';
    fprintf(stderr, "MPI-IO Error: Could not open file \"%s\": %s\n",
            fn.c_str(), buf);
    fatal_error();
  }
  ret = MPI_File_set_view(f, pref * sizeof(T), MPI_T, MPI_T,
                          const_cast<char *>("native"), MPI_INFO_NULL);

  ret |= MPI_File_read_all(f, arr, len, MPI_T, MPI_STATUS_IGNORE);
  MPI_File_close(&f);
  if (ret) {
    fprintf(stderr, "MPI-IO Error: Could not read file \"%s\".\n", fn.c_str());
    fatal_error();
  }
}

/**
 * @brief Read the header file and return the first value.
 * To be called by all processes.
 *
 * @param fn Filename of the head file
 * @param rank The rank of the current process in @c MPI_COMM_WORLD
 */
static unsigned read_head(const std::string &fn, int rank) {
  unsigned n_fields = 0u;
  FILE *f = nullptr;
  if (rank == 0) {
    if (!(f = fopen(fn.c_str(), "rb"))) {
      fprintf(stderr, "MPI-IO: Could not open \"%s\".\n", fn.c_str());
      fatal_error();
    }
    if (fread(static_cast<void *>(&n_fields), sizeof(unsigned), 1, f) != 1) {
      fprintf(stderr, "MPI-IO: Read on \"%s\" failed.\n", fn.c_str());
      fatal_error();
    }
  }
  MPI_Bcast(&n_fields, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  if (f) {
    fclose(f);
  }
  return n_fields;
}

/**
 * @brief Read the pref file.
 * Needs to be called by all processes.
 *
 * @param fn The file name of the prefs file
 * @param rank The rank of the current process in @c MPI_COMM_WORLD
 * @param size The size of @c MPI_COMM_WORLD
 * @param nglobalpart The global amount of particles
 * @return The prefix and the local number of particles.
 */
static std::tuple<int, int> read_prefs(const std::string &fn, int rank,
                                       int size, int nglobalpart) {
  int pref = 0;
  int nlocalpart = 0;
  mpiio_read_array<int>(fn, &pref, 1, rank, MPI_INT);
  if (rank > 0)
    MPI_Send(&pref, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD);
  if (rank < size - 1)
    MPI_Recv(&nlocalpart, 1, MPI_INT, rank + 1, MPI_ANY_TAG, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
  else
    nlocalpart = nglobalpart;
  nlocalpart -= pref;
  return {pref, nlocalpart};
}

void mpi_mpiio_common_read(const std::string &prefix, unsigned fields) {
  cell_structure.remove_all_particles();

  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  auto const nproc = get_num_elem(prefix + ".pref", sizeof(int));
  auto const nglobalpart = get_num_elem(prefix + ".id", sizeof(int));

  if (rank == 0 && nproc != size) {
    fprintf(stderr, "MPI-IO Error: Trying to read a file with a different COMM "
                    "size than at point of writing.\n");
    fatal_error();
  }

  // 1.head on head node:
  // Read head to determine fields at time of writing.
  // Compare this var to the current fields.
  auto const avail_fields = read_head(prefix + ".head", rank);
  if (rank == 0 && (fields & avail_fields) != fields) {
    fprintf(stderr,
            "MPI-IO Error: Requesting to read fields which were not dumped.\n");
    fatal_error();
  }

  // 1.pref on all nodes:
  // Read own prefix (1 int at prefix rank).
  // Communicate own prefix to rank-1
  // Determine nlocalpart (prefix of rank+1 - own prefix) on every node.
  int pref, nlocalpart;
  std::tie(pref, nlocalpart) =
      read_prefs(prefix + ".pref", rank, size, nglobalpart);

  std::vector<Particle> particles(nlocalpart);

  {
    // 1.id on all nodes:
    // Read nlocalpart ints at defined prefix.
    std::vector<int> id(nlocalpart);
    mpiio_read_array<int>(prefix + ".id", id.data(), nlocalpart, pref, MPI_INT);

    for (int i = 0; i < nlocalpart; ++i) {
      particles[i].p.identity = id[i];
    }
  }

  if (fields & MPIIO_OUT_POS) {
    // 1.pos on all nodes:
    // Read nlocalpart * 3 doubles at defined prefix * 3
    std::vector<double> pos(3 * nlocalpart);
    mpiio_read_array<double>(prefix + ".pos", pos.data(), 3 * nlocalpart,
                             3 * pref, MPI_DOUBLE);

    for (int i = 0; i < nlocalpart; ++i) {
      particles[i].r.p = {pos[3 * i + 0], pos[3 * i + 1], pos[3 * i + 2]};
    }
  }

  if (fields & MPIIO_OUT_TYP) {
    // 1.type on all nodes:
    // Read nlocalpart ints at defined prefix.
    std::vector<int> type(nlocalpart);
    mpiio_read_array<int>(prefix + ".type", type.data(), nlocalpart, pref,
                          MPI_INT);

    for (int i = 0; i < nlocalpart; ++i)
      particles[i].p.type = type[i];
  }

  if (fields & MPIIO_OUT_VEL) {
    // 1.vel on all nodes:
    // Read nlocalpart * 3 doubles at defined prefix * 3
    std::vector<double> vel(3 * nlocalpart);
    mpiio_read_array<double>(prefix + ".vel", vel.data(), 3 * nlocalpart,
                             3 * pref, MPI_DOUBLE);

    for (int i = 0; i < nlocalpart; ++i)
      particles[i].m.v = {vel[3 * i + 0], vel[3 * i + 1], vel[3 * i + 2]};
  }

  if (fields & MPIIO_OUT_BND) {
    // 1.boff
    // 1 int per process
    int bonds_size = 0;
    mpiio_read_array<int>(prefix + ".boff", &bonds_size, 1, rank, MPI_INT);
    int bpref = 0;
    MPI_Exscan(&bonds_size, &bpref, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    // 1.bond
    // nlocalbonds ints per process
    std::vector<char> bond(bonds_size);
    mpiio_read_array<char>(prefix + ".bond", bond.data(), bonds_size, bpref,
                           MPI_CHAR);

    namespace io = boost::iostreams;
    io::array_source src(bond.data(), bond.size());
    io::stream<io::array_source> ss(src);
    boost::archive::binary_iarchive ia(ss);

    for (auto &p : particles) {
      ia >> p.bonds();
    }
  }

  for (auto &p : particles) {
    cell_structure.add_particle(std::move(p));
  }
}
} // namespace Mpiio
