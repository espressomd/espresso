/*
 * Copyright (C) 2010-2022 The ESPResSo project
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
 *   <tt>rank0 --- rank1 --- rank2 ...</tt>
 *   where each rank dumps its scalars in the ordering of the particles.
 * - Vector arrays are written in the rank ordering like scalar arrays.
 *   The ordering of the vector data is: <tt>v[0] v[1] v[2]</tt>, so the data
 *   looks like this:
 *   <tt>v1[0] v1[1] v1[2] v2[0] v2[1] v2[2] v3[0] ...</tt>
 *
 * To be able to determine the rank boundaries (a multiple of
 * @c nlocalparts), the file 1.pref is written, which dumps the partial
 * sum of @c nlocalparts, i.e. the prefixes in scalar arrays:
 * - 1.prefs looks like this:
 *   <tt>0 nlocalpats_rank0 nlocalparts_rank0+nlocalparts_rank1 ...</tt>
 *
 * Bonds are dumped as two arrays, namely 1.bond which stores the
 * bonding partners of the particles and 1.boff which stores the
 * iteration indices for each particle.
 * - 1.boff is a scalar array of size <tt>(nlocalpart + 1)</tt> per rank.
 * - The last element (at index @c nlocalpart) of 1.boff's subpart
 *   <tt>[rank * (nlocalpart + 1) : (rank + 1) * (nlocalpart + 1)]</tt>
 *   determines the number of bonds for processor @c rank.
 * - In this subarray one can find the bonding partners of particle
 *   <tt>id[i]</tt>. The iteration indices for local part of 1.bonds are:
 *   <tt>subarray[i] : subarray[i+1]</tt>
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

#include <algorithm>
#include <cassert>
#include <cerrno>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <sstream>
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
 * @param msg          Custom error message
 * @param fn           File path
 * @param extra        Extra context
 */
static bool fatal_error(char const *msg, std::string const &fn = "",
                        std::string const &extra = "") {
  std::stringstream what;
  what << "MPI-IO Error: " << msg;
  if (not fn.empty()) {
    what << " \"" << fn << "\"";
  }
  if (not extra.empty()) {
    what << " :" << extra;
  }
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size == 1) {
    throw std::runtime_error(what.str());
  }
  fprintf(stderr, "%s\n", what.str().c_str());
  errexit();
  return false;
}

/**
 * @brief Fatal error handler that closes an open file and queries the
 * message associated with an MPI error code.
 * On 1 MPI rank the error is recoverable and an exception is thrown.
 * On more than 1 MPI rank the error is not recoverable.
 * @param msg          Custom error message
 * @param fn           File path
 * @param fp           File handle
 * @param errnum       MPI error code
 */
static bool fatal_error(char const *msg, std::string const &fn, MPI_File *fp,
                        int errnum) {
  // get MPI error message
  char buf[MPI_MAX_ERROR_STRING];
  int buf_len;
  MPI_Error_string(errnum, buf, &buf_len);
  buf[buf_len] = '\0';
  // close file handle
  if (fp) {
    MPI_File_close(fp);
  }
  return fatal_error(msg, fn, buf);
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
static void mpiio_dump_array(const std::string &fn, T const *arr,
                             std::size_t len, std::size_t pref,
                             MPI_Datatype MPI_T) {
  MPI_File f;
  int ret;
  ret = MPI_File_open(MPI_COMM_WORLD, const_cast<char *>(fn.c_str()),
                      // MPI_MODE_EXCL: Prohibit overwriting
                      MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_EXCL,
                      MPI_INFO_NULL, &f);
  if (ret) {
    fatal_error("Could not open file", fn, &f, ret);
  }
  auto const offset =
      static_cast<MPI_Offset>(pref) * static_cast<MPI_Offset>(sizeof(T));
  ret = MPI_File_set_view(f, offset, MPI_T, MPI_T, const_cast<char *>("native"),
                          MPI_INFO_NULL);
  ret |= MPI_File_write_all(f, arr, static_cast<int>(len), MPI_T,
                            MPI_STATUS_IGNORE);
  static_cast<void>(ret and fatal_error("Could not write file", fn, &f, ret));
  MPI_File_close(&f);
}

/**
 * @brief Calculate the file offset on the local node.
 * @param n_items  Number of items on the local node.
 * @return The number of items on all nodes with lower rank.
 */
static unsigned long mpi_calculate_file_offset(unsigned long n_items) {
  unsigned long offset = 0ul;
  MPI_Exscan(&n_items, &offset, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  return offset;
}

/**
 * @brief Dump the fields and bond information.
 * To be called by the head node only.
 *
 * @param fn The filename to write to
 * @param fields The dumped fields
 */
static void dump_info(const std::string &fn, unsigned fields) {
  // MPI-IO requires consecutive bond ids
  auto const nbonds = bonded_ia_params.size();
  assert(static_cast<std::size_t>(bonded_ia_params.get_next_key()) == nbonds);

  FILE *f = fopen(fn.c_str(), "wb");
  if (!f) {
    fatal_error("Could not open file", fn);
  }
  static std::vector<int> npartners;
  bool success = (fwrite(&fields, sizeof(fields), 1u, f) == 1);
  // Pack the necessary information of bonded_ia_params:
  // The number of partners. This is needed to interpret the bond IntList.
  if (nbonds > npartners.size())
    npartners.resize(nbonds);

  auto npartners_it = npartners.begin();
  for (int i = 0; i < bonded_ia_params.get_next_key(); ++i, ++npartners_it) {
    *npartners_it = number_of_partners(*bonded_ia_params.at(i));
  }
  success = success && (fwrite(&nbonds, sizeof(std::size_t), 1u, f) == 1);
  success =
      success && (fwrite(npartners.data(), sizeof(int), nbonds, f) == nbonds);
  fclose(f);
  static_cast<void>(success or fatal_error("Could not write file", fn));
}

void mpi_mpiio_common_write(const std::string &prefix, unsigned fields,
                            const ParticleRange &particles) {
  auto const nlocalpart = static_cast<unsigned long>(particles.size());
  auto const offset = mpi_calculate_file_offset(nlocalpart);
  // Keep static buffers in order to avoid allocating them on every
  // function call
  static std::vector<double> pos, vel;
  static std::vector<int> id, type;

  // Realloc static buffers if necessary
  if (nlocalpart > id.size())
    id.resize(nlocalpart);
  if (fields & MPIIO_OUT_POS && 3ul * nlocalpart > pos.size())
    pos.resize(3ul * nlocalpart);
  if (fields & MPIIO_OUT_VEL && 3ul * nlocalpart > vel.size())
    vel.resize(3ul * nlocalpart);
  if (fields & MPIIO_OUT_TYP && nlocalpart > type.size())
    type.resize(nlocalpart);

  // Pack the necessary information
  auto id_it = id.begin();
  auto type_it = type.begin();
  auto pos_it = pos.begin();
  auto vel_it = vel.begin();
  for (auto const &p : particles) {
    *id_it = p.id();
    ++id_it;
    if (fields & MPIIO_OUT_POS) {
      std::copy_n(std::begin(p.pos()), 3u, pos_it);
      pos_it += 3u;
    }
    if (fields & MPIIO_OUT_VEL) {
      std::copy_n(std::begin(p.v()), 3u, vel_it);
      vel_it += 3u;
    }
    if (fields & MPIIO_OUT_TYP) {
      *type_it = p.type();
      ++type_it;
    }
  }

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0)
    dump_info(prefix + ".head", fields);
  auto const pref_offset = static_cast<unsigned long>(rank);
  mpiio_dump_array<unsigned long>(prefix + ".pref", &offset, 1ul, pref_offset,
                                  MPI_UNSIGNED_LONG);
  mpiio_dump_array<int>(prefix + ".id", id.data(), nlocalpart, offset, MPI_INT);
  if (fields & MPIIO_OUT_POS)
    mpiio_dump_array<double>(prefix + ".pos", pos.data(), 3ul * nlocalpart,
                             3ul * offset, MPI_DOUBLE);
  if (fields & MPIIO_OUT_VEL)
    mpiio_dump_array<double>(prefix + ".vel", vel.data(), 3ul * nlocalpart,
                             3ul * offset, MPI_DOUBLE);
  if (fields & MPIIO_OUT_TYP)
    mpiio_dump_array<int>(prefix + ".type", type.data(), nlocalpart, offset,
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
    auto const bonds_size = static_cast<unsigned long>(bonds.size());
    auto const bonds_offset = mpi_calculate_file_offset(bonds_size);

    mpiio_dump_array<unsigned long>(prefix + ".boff", &bonds_size, 1ul,
                                    pref_offset, MPI_UNSIGNED_LONG);
    mpiio_dump_array<char>(prefix + ".bond", bonds.data(), bonds.size(),
                           bonds_offset, MPI_CHAR);
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
static unsigned long get_num_elem(const std::string &fn, std::size_t elem_sz) {
  // Could also be done via MPI_File_open, MPI_File_get_size,
  // MPI_File_close.
  struct stat st;
  errno = 0;
  if (stat(fn.c_str(), &st) != 0) {
    auto const reason = strerror(errno);
    fatal_error("Could not get file size of", fn, reason);
  }
  return static_cast<unsigned long>(st.st_size) / elem_sz;
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
    fatal_error("Could not open file", fn, &f, ret);
  }
  auto const offset =
      static_cast<MPI_Offset>(pref) * static_cast<MPI_Offset>(sizeof(T));
  ret = MPI_File_set_view(f, offset, MPI_T, MPI_T, const_cast<char *>("native"),
                          MPI_INFO_NULL);

  ret |= MPI_File_read_all(f, arr, static_cast<int>(len), MPI_T,
                           MPI_STATUS_IGNORE);
  static_cast<void>(ret and fatal_error("Could not read file", fn, &f, ret));
  MPI_File_close(&f);
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
    f = fopen(fn.c_str(), "rb");
    static_cast<void>(not f and fatal_error("Could not open file", fn));
    auto const n = fread(static_cast<void *>(&n_fields), sizeof n_fields, 1, f);
    static_cast<void>((n == 1) or fatal_error("Could not read file", fn));
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
static std::tuple<unsigned long, unsigned long>
read_prefs(const std::string &fn, int rank, int size,
           unsigned long nglobalpart) {
  auto const pref_offset = static_cast<unsigned long>(rank);
  unsigned long pref = 0ul;
  unsigned long nlocalpart = 0ul;
  mpiio_read_array<unsigned long>(fn, &pref, 1ul, pref_offset,
                                  MPI_UNSIGNED_LONG);
  if (rank > 0)
    MPI_Send(&pref, 1, MPI_UNSIGNED_LONG, rank - 1, 0, MPI_COMM_WORLD);
  if (rank < size - 1)
    MPI_Recv(&nlocalpart, 1, MPI_UNSIGNED_LONG, rank + 1, MPI_ANY_TAG,
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
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
  auto const nproc = get_num_elem(prefix + ".pref", sizeof(unsigned long));
  auto const nglobalpart = get_num_elem(prefix + ".id", sizeof(int));

  if (rank == 0 && nproc != static_cast<unsigned long>(size)) {
    fatal_error("Trying to read a file with a different COMM "
                "size than at point of writing.");
  }

  // 1.head on head node:
  // Read head to determine fields at time of writing.
  // Compare this var to the current fields.
  auto const avail_fields = read_head(prefix + ".head", rank);
  if (rank == 0 && (fields & avail_fields) != fields) {
    fatal_error("Requesting to read fields which were not dumped.");
  }

  // 1.pref on all nodes:
  // Read own prefix (1 int at prefix rank).
  // Communicate own prefix to rank-1
  // Determine nlocalpart (prefix of rank+1 - own prefix) on every node.
  auto const [pref, nlocalpart] =
      read_prefs(prefix + ".pref", rank, size, nglobalpart);

  std::vector<Particle> particles(nlocalpart);

  {
    // 1.id on all nodes:
    // Read nlocalpart ints at defined prefix.
    std::vector<int> id(nlocalpart);
    auto id_it = id.begin();
    mpiio_read_array<int>(prefix + ".id", id.data(), nlocalpart, pref, MPI_INT);

    for (auto &p : particles) {
      p.id() = *id_it;
      ++id_it;
    }
  }

  if (fields & MPIIO_OUT_POS) {
    // 1.pos on all nodes:
    // Read nlocalpart * 3 doubles at defined prefix * 3
    std::vector<double> pos(3ul * nlocalpart);
    auto pos_it = pos.begin();
    mpiio_read_array<double>(prefix + ".pos", pos.data(), 3ul * nlocalpart,
                             3ul * pref, MPI_DOUBLE);

    for (auto &p : particles) {
      std::copy_n(pos_it, 3u, std::begin(p.pos()));
      pos_it += 3u;
    }
  }

  if (fields & MPIIO_OUT_TYP) {
    // 1.type on all nodes:
    // Read nlocalpart ints at defined prefix.
    std::vector<int> type(nlocalpart);
    auto type_it = type.begin();
    mpiio_read_array<int>(prefix + ".type", type.data(), nlocalpart, pref,
                          MPI_INT);

    for (auto &p : particles) {
      p.type() = *type_it;
      ++type_it;
    }
  }

  if (fields & MPIIO_OUT_VEL) {
    // 1.vel on all nodes:
    // Read nlocalpart * 3 doubles at defined prefix * 3
    std::vector<double> vel(3ul * nlocalpart);
    auto vel_it = vel.begin();
    mpiio_read_array<double>(prefix + ".vel", vel.data(), 3ul * nlocalpart,
                             3ul * pref, MPI_DOUBLE);

    for (auto &p : particles) {
      std::copy_n(vel_it, 3u, std::begin(p.v()));
      vel_it += 3u;
    }
  }

  if (fields & MPIIO_OUT_BND) {
    // 1.boff
    // 1 long int per process
    auto const pref_offset = static_cast<unsigned long>(rank);
    unsigned long bonds_size = 0u;
    mpiio_read_array<unsigned long>(prefix + ".boff", &bonds_size, 1ul,
                                    pref_offset, MPI_UNSIGNED_LONG);
    auto const bonds_offset = mpi_calculate_file_offset(bonds_size);

    // 1.bond
    // nlocalbonds ints per process
    std::vector<char> bond(bonds_size);
    mpiio_read_array<char>(prefix + ".bond", bond.data(), bonds_size,
                           bonds_offset, MPI_CHAR);

    boost::iostreams::array_source src(bond.data(), bond.size());
    boost::iostreams::stream<boost::iostreams::array_source> ss(src);
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
