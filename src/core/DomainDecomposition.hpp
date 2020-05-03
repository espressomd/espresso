#ifndef ESPRESSO_DOMAIN_DECOMPOSITION_HPP
#define ESPRESSO_DOMAIN_DECOMPOSITION_HPP

#include "ParticleDecomposition.hpp"

#include "BoxGeometry.hpp"
#include "LocalBox.hpp"

#include <boost/range/numeric.hpp>

#include <utils/index.hpp>
#include <utils/mpi/cart_comm.hpp>

/** Structure containing the information about the cell grid used for domain
 *  decomposition.
 */
struct DomainDecomposition : public ParticleDecomposition {
  /** Grind dimensions per node. */
  Utils::Vector3i cell_grid = {};
  /** cell size. */
  Utils::Vector3d cell_size = {};

private:
  /** Offset in global grid */
  Utils::Vector3i cell_offset = {};
  /** linked cell grid with ghost frame. */
  Utils::Vector3i ghost_cell_grid = {};
  /** inverse cell size = \see DomainDecomposition::cell_size ^ -1. */
  Utils::Vector3d inv_cell_size = {};

  boost::mpi::communicator comm;
  BoxGeometry box_geo;
  LocalBox<double> local_geo;
  std::vector<Cell> cells;
  std::vector<Cell *> m_local_cells;
  std::vector<Cell *> m_ghost_cells;
  GhostCommunicator m_exchange_ghosts_comm;
  GhostCommunicator m_collect_ghost_force_comm;

public:
  DomainDecomposition() = default;
  DomainDecomposition(const boost::mpi::communicator &comm, double range,
                      const BoxGeometry &box_geo,
                      const LocalBox<double> &local_geo);

  GhostCommunicator const &exchange_ghosts_comm() const override {
    return m_exchange_ghosts_comm;
  }
  GhostCommunicator const &collect_ghost_force_comm() const override {
    return m_collect_ghost_force_comm;
  };

  Utils::Span<Cell *> local_cells() override {
    return Utils::make_span(m_local_cells);
  }
  Utils::Span<Cell *> ghost_cells() override {
    return Utils::make_span(m_ghost_cells);
  }

  Cell *particle_to_cell(Particle const &p) override {
    return position_to_cell(p.r.p);
  }

  bool minimum_image_distance() const override { return false; }

  void resort(bool global, ParticleList &pl,
              std::vector<Cell *> &modified_cells) override;

private:
  /** Fill local_cells list and ghost_cells list for use with domain
   *  decomposition.  \ref cells::cells is assumed to be a 3d grid with size
   *  \ref DomainDecomposition::ghost_cell_grid.
   */
  void mark_cells();

  /** Fill a communication cell pointer list. Fill the cell pointers of
   *  all cells which are inside a rectangular subgrid of the 3D cell
   *  grid (\ref DomainDecomposition::ghost_cell_grid) starting from the
   *  lower left corner lc up to the high top corner hc. The cell
   *  pointer list part_lists must already be large enough.
   *  \param part_lists  List of cell pointers to store the result.
   *  \param lc          lower left corner of the subgrid.
   *  \param hc          high up corner of the subgrid.
   */
  void fill_comm_cell_lists(ParticleList **part_lists,
                            Utils::Vector3i const &lc,
                            Utils::Vector3i const &hc);

  /** @brief Maximal interaction range supported with
   *         the current geometry and node grid.
   * @return Per-direction maximal range.
   */
  Utils::Vector3d max_range() const override;

  int calc_processor_min_num_cells() const;

private:
  Cell *position_to_cell(const Utils::Vector3d &pos);

  /**
   * @brief Move particles into the cell system if it belongs to this node.
   *
   * Moves all particles from src into the local cell
   * system if they do belong here. Otherwise the
   * particles are moved into rest.
   *
   * @param src Particles to move.
   * @param rest Output list for left-over particles.
   * @param modified_cells Local cells that were touched.
   */
  void move_if_local(ParticleList &src, ParticleList &rest,
                     std::vector<Cell *> &modified_cells);

  /**
   * @brief Split particle list by direction.
   *
   * Moves all particles from src into left
   * and right depending if they belong to
   * the left or right side from local node
   * in direction dir.
   *
   * @param src Particles to sort.
   * @param left Particles that should go to the left
   * @param right Particles that should go to the right
   * @param dir Direction to consider.
   */
  void move_left_or_right(ParticleList &src, ParticleList &left,
                          ParticleList &right, int dir) const;

  /**
   * @brief One round of particle exchange with the next neighbors.
   *
   * @param[in] pl Particle on the move
   * @param[out] modified_cells Cells that got touched.
   */
  void exchange_neighbors(ParticleList &pl,
                          std::vector<Cell *> &modified_cells);

  /**
   *  @brief Calculate cell grid dimensions, cell sizes and number of cells.
   *
   *  Calculates the cell grid, based on \ref local_geo and \p range.
   *  If the number of cells is larger than \ref max_num_cells,
   *  it increases max_range until the number of cells is
   *  smaller or equal \ref max_num_cells. It sets:
   *  \ref DomainDecomposition::cell_grid,
   *  \ref DomainDecomposition::ghost_cell_grid,
   *  \ref DomainDecomposition::cell_size, and
   *  \ref DomainDecomposition::inv_cell_size.
   *
   *  @param range Required interacting range. All pairs closer
   *         than this distance are found.
   */
  void create_cell_grid(double range);

  /** Init cell interactions for cell system domain decomposition.
   * initializes the interacting neighbor cell list of a cell The
   * created list of interacting neighbor cells is used by the Verlet
   * algorithm.
   */
  void init_cell_interactions();

  /** Create communicators for cell structure domain decomposition. (see \ref
   *  GhostCommunicator)
   */
  GhostCommunicator prepare_comm();

  /** update the 'shift' member of those GhostCommunicators, which use
   *  that value to speed up the folding process of its ghost members
   *  (see \ref dd_prepare_comm for the original).
   */
  void update_communicators_w_boxl();

  bool on_geometry_change(bool fast, double range, const BoxGeometry &box_geo,
                          const LocalBox<double> &local_geo);

public:
  /** Maximal number of cells per node. In order to avoid memory
   *  problems due to the cell grid one has to specify the maximal
   *  number of cells. If the number of cells is larger
   *  than max_num_cells the cell grid is reduced.
   *  max_num_cells has to be larger than 27, e.g. one inner cell.
   */
  static constexpr int max_num_cells = 32768;
};

#endif
