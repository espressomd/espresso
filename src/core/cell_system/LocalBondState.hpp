#pragma once
#include <config/config.hpp>

#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM

#include <Kokkos_Core.hpp>

struct LocalBondState {
  using PairBondlistType = Kokkos::View<int *[2], Kokkos::LayoutRight>;
  using PairBondIDType = Kokkos::View<int *, Kokkos::LayoutRight>;
  using AngleBondlistType = Kokkos::View<int *[3], Kokkos::LayoutRight>;
  using AngleBondIDType = Kokkos::View<int *, Kokkos::LayoutRight>;
  using DihedralBondlistType = Kokkos::View<int *[4], Kokkos::LayoutRight>;
  using DihedralBondIDType = Kokkos::View<int *, Kokkos::LayoutRight>;

  // Bond counts
  int pair_count = 0;
  int angle_count = 0;
  int dihedral_count = 0;

  // Kokkos bond lists
  PairBondlistType pair_list;
  PairBondIDType pair_ids;
  AngleBondlistType angle_list;
  AngleBondIDType angle_ids;
  DihedralBondlistType dihedral_list;
  DihedralBondIDType dihedral_ids;

#ifdef ESPRESSO_COLLISION_DETECTION
  std::vector<int> new_pair_list, new_pair_ids;
  std::vector<int> new_angle_list, new_angle_ids;
  std::vector<int> new_dihedral_list, new_dihedral_ids;
#endif

  void reset_counts() {
    pair_count = 0;
    angle_count = 0;
    dihedral_count = 0;
  }

  void set_counts(int p, int a, int d) {
    pair_count = p;
    angle_count = a;
    dihedral_count = d;
  }

  /// Allocate or reallocate all Kokkos Views to current counts.
  void allocate();

  /// Full clear — deallocates Views
  void clear();

  /// Light reset — only reset counts + collision vectors
  void reset();

#ifdef ESPRESSO_COLLISION_DETECTION
  void clear_new_bonds();
  void add_new_bond(int bond_id, std::vector<int> const &particle_ids,
                    Kokkos::View<int *> const &id_to_index);
  void rebuild();
#endif
};

#endif // ESPRESSO_SHARED_MEMORY_PARALLELISM
