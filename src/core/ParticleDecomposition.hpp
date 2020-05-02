#ifndef ESPRESSO_PARTICLE_DECOMPOSITION_HPP
#define ESPRESSO_PARTICLE_DECOMPOSITION_HPP

#include "BoxGeometry.hpp"
#include "Cell.hpp"
#include "LocalBox.hpp"
#include "ghosts.hpp"

#include <utils/Span.hpp>
#include <utils/Vector.hpp>

class ParticleDecomposition {
public:
  /**
   * @brief Resort particles.
   *
   * @param global_flag[in] Expect particles to be displaced by more than a
   * local box size.
   * @param displaced_parts[inout] Particles to be sorted
   * @param modified_cells[out] Cells that have been touched.
   */
  virtual void resort(bool global_flag, ParticleList &displaced_parts,
                      std::vector<Cell *> &modified_cells) = 0;

  /**
   * @brief Communicator for updateing ghosts from the real particles.
   */
  virtual GhostCommunicator const &exchange_ghosts_comm() const = 0;
  /**
   * @brief Communicator for force reduction.
   */
  virtual GhostCommunicator const &collect_ghost_force_comm() const = 0;

  /**
   * @brief Get pointer to local cells.
   *
   * Local cells are cells that contain particles
   * that are owned by this node.
   *
   * @return List of local cells.
   */
  virtual Utils::Span<Cell *> local_cells() = 0;

  /**
   * @brief Get pointer to local cells.
   *
   * Ghost cells are cells that contain particles
   * that are owned by different nodes but interact
   * with particles on this node.
   *
   * @return List of ghost cells.
   */
  virtual Utils::Span<Cell *> ghost_cells() = 0;

  /**
   * @brief Determine which cell a particle id belongs to.
   *
   * @param p Particle to find cell for.
   * @return Pointer to cell or nullptr if not local.
   */
  virtual Cell *particle_to_cell(Particle const &p) = 0;

  /**
   * @brief Maximum supported cutoff.
   */
  virtual Utils::Vector3d max_range() const = 0;
  /**
   * @brief Return true if minimum image convention is
   *        needed for distance calculation. */
  virtual bool minimum_image_distance() const = 0;

  /**
   * @brief Try to change the box geometry.
   *
   * This method can be used to change the geometry properties
   * of the decomposition, e.g. the global and local box size
   * as well as the interaction range. This call can fail, in
   * which case false is returned, and the decomposition remains
   * unchanged. In this case all pointers into the decomposition
   * stay valid.
   *
   *  @param fast If true return asap.
   *  @param range Desired interaction range
   *  @param box_geo New box geometry.
   *  @param local_geo New local box.
   *  @return If the change was possible.
   */
  virtual bool on_geometry_change(bool fast, double range,
                                  const BoxGeometry &box_geo,
                                  const LocalBox<double> &local_geo) = 0;

  virtual ~ParticleDecomposition() = default;
};

#endif
