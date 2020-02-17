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
#ifndef _GHOSTS_H
#define _GHOSTS_H
/** \file
 *  Ghost particles and particle exchange.
 *
 *  In this file you find everything concerning the exchange of
 *  particle data (particles, ghosts, positions and forces) for short
 *  range interactions between the spacial domains of neighbouring
 *  nodes.
 *
 *  <h2> How does this work </h2>
 *  The ghost communication transfers data from cells on one node to cells on
 *  another node during the integration process. Note that data can only be
 *  transferred from one cell to one other, and the contents of the other cell
 *  will be overwritten. This communication is invoked by the integrator, at
 *  four different times in an integration step. These steps are reflected in
 *  @ref cell_structure structure, since they are treated differently by
 *  different cell systems:
 *  - @ref CellStructure::m_ghost_cells are used to transfer the cell sizes,
 *    i.e., make sure that for all later transfers the target cell has the same
 *    size as the source cell.
 *  - @ref CellStructure::exchange_ghosts_comm sets up all information on the
 *    ghosts that is necessary. Normally transfers the (shifted) position and
 *    particle properties.
 *  - @ref cells_update_ghosts is used to update the particle properties if no
 *    particles have been moved between cells. Therefore only the positions are
 *    transferred, otherwise this normally looks pretty much like the
 *    @ref CellStructure::exchange_ghosts_comm.
 *  - @ref CellStructure::collect_ghost_force_comm finally is used to transfer
 *    back forces that were exerted on a ghost particle. They are simply added
 *    again to the force of the real particle. The communication process is
 *    therefore inverted with respect to
 *    @ref CellStructure::exchange_ghosts_comm and @ref cells_update_ghosts,
 *    i.e., sending is replaced by receiving and the other way round.
 *
 *  The particle data that has to be transferred, and especially from where to
 *  where, heavily depends on the cell system. In ESPResSo, this is abstracted
 *  in form of ghost communicators (\ref GhostCommunicator) and ghost
 *  communications (\ref GhostCommunication). The ghost communicators represent
 *  the four communications above and consist of the data to transfer (which is
 *  determined by their type) and a list of ghost communications. The data
 *  types are described by the particle data classes:
 *  - @ref GHOSTTRANS_PROPRTS transfers the @ref ParticleProperties
 *  - @ref GHOSTTRANS_POSITION transfers the @ref ParticlePosition
 *  - @ref GHOSTTRANS_MOMENTUM transfers the @ref ParticleMomentum
 *  - @ref GHOSTTRANS_FORCE transfers the @ref ParticleForce
 *  - @ref GHOSTTRANS_PARTNUM transfers the cell sizes
 *
 *  Each ghost communication describes a single communication of the local with
 *  another node (or all other nodes). The data transferred can be any number
 *  of cells, there are five communication types:
 *  - @ref GHOST_SEND sends data to one other node
 *  - @ref GHOST_RECV recvs data from one other node. In the case of
 *    forces, they are added up, everything else is overwritten
 *  - @ref GHOST_BCST sends data to all nodes
 *  - @ref GHOST_RDCE recvs data from all nodes. In the case of
 *    forces, they are added up, everything else is overwritten
 *  - @ref GHOST_LOCL transfer data from a local cell to another local cell.
 *    In this case, the first half of the cells are the sending cells, the
 *    other half are the corresponding receivers.
 *
 *  Note that for the first four communications you have to make sure that
 *  when one node sends to another, that the sender has @ref GHOST_SEND
 *  and the receiver @ref GHOST_RECV. In the case of @ref GHOST_BCST resp.
 *  @ref GHOST_RDCE, all nodes have to have the same communication type
 *  and the same master sender/receiver (just like the MPI commands).
 *
 *  The ghost communicators are created in the init routines of the cell
 *  systems, therefore have a look at @ref dd_topology_init or
 *  @ref nsq_topology_init for further details.
 */
#include "ParticleList.hpp"

#include <boost/mpi/communicator.hpp>

/** \name Transfer types, for \ref GhostCommunicator::type */
/************************************************************/
/*@{*/

enum : unsigned {
  /// send to a single node
  GHOST_SEND = 0,
  /// recv from a single node
  GHOST_RECV = 1,
  /// broadcast, the node entry gives the sender
  GHOST_BCST = 2,
  /// reduce, the node entry gives the receiver
  GHOST_RDCE = 3,
  /// transfer data from cell to cell on this node
  GHOST_LOCL = 4
};

/*@}*/

/** Transfer data classes, for \ref ghost_communicator */
enum : unsigned {
  GHOSTTRANS_NONE = 0u,
  /// transfer \ref ParticleProperties
  GHOSTTRANS_PROPRTS = 1u,
  /// transfer \ref ParticlePosition
  GHOSTTRANS_POSITION = 2u,
  /// transfer \ref ParticleMomentum
  GHOSTTRANS_MOMENTUM = 8u,
  /// transfer \ref ParticleForce
  GHOSTTRANS_FORCE = 16u,
  /// resize the receiver particle arrays to the size of the senders
  GHOSTTRANS_PARTNUM = 64u,
  GHOSTTRANS_BONDS = 128u
};

/** \name Data Types */
/************************************************************/
/*@{*/

/**
 * Class that stores marshalled data for ghost communications.
 * To store and retrieve data, use the adapter classes below.
 */
class CommBuf {
public:
  /** Returns a pointer to the non-bond storage.
   */
  char *data() { return buf.data(); }

  /** Returns the number of elements in the non-bond storage.
   */
  size_t size() { return buf.size(); }

  /** Resizes the underlying storage s.t. the object is capable
   * of holding "new_size" chars.
   * @param new_size new size
   */
  void resize(size_t new_size) { buf.resize(new_size); }

  /** Returns a reference to the bond storage.
   */
  std::vector<int> &bonds() { return bondbuf; }

private:
  std::vector<char> buf;    //< Buffer for everything but bonds
  std::vector<int> bondbuf; //< Buffer for bond lists
};

struct GhostCommunication {
  unsigned type = GHOSTTRANS_NONE;
  boost::mpi::communicator comm;
  /** Node to communicate with (to use with all MPI operations). */
  int send_to = -1;
  int recv_from = -1;

  /** Pointer array to particle lists to communicate. */
  std::vector<ParticleList *> send_lists = {};
  std::vector<ParticleList *> recv_lists = {};

  /** Position shift for ghost particles. The shift is done on the sender side.
   */
  Utils::Vector3d shift = {};

  /** Buffers for MPI serialization */
  mutable CommBuf send_buffer, recv_buffer;
};

/** Properties for a ghost communication. A ghost communication is defined */
struct GhostCommunicator {
  GhostCommunicator() = default;
  explicit GhostCommunicator(std::vector<GhostCommunication> communications)
      : m_communications(std::move(communications)) {}

  /** List of ghost communications. */
  const std::vector<GhostCommunication> &communications() const & {
    return m_communications;
  }
  std::vector<GhostCommunication> &&communications() && {
    return std::move(m_communications);
  }

private:
  std::vector<GhostCommunication> m_communications;
};

/*@}*/

/** \name Exported Functions */
/************************************************************/
/*@{*/

/**
 * @brief Do a ghost communication with caller specified data parts.
 */
void ghost_communicator(const GhostCommunicator *gcr, unsigned int data_parts);

/*@}*/

#endif
