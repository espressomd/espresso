#
# Copyright (C) 2025 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

import espressomd
import numpy as np
import argparse
import time


def benchmark_operation(operation_func, n_iter):
    """Execute an operation multiple times and return timing statistics."""
    timings = []
    for _ in range(n_iter):
        tick = time.time()
        operation_func()
        tock = time.time()
        timings.append(tock - tick)
    return np.mean(timings), np.std(timings)


def print_section_header(title):
    """Print a formatted section header."""
    print("=" * 60)
    print(title)
    print("=" * 60)


def filter_available_properties(properties, particle):
    """Filter properties to only those available for the given particle."""
    return [prop for prop in properties if hasattr(particle, prop)]


parser = argparse.ArgumentParser(
    description="Benchmark particle creation and property access.")
parser.add_argument("--particles_per_core", metavar="N", action="store",
                    type=int, default=1000, required=False,
                    help="Number of particles per core (default: 1000)")
parser.add_argument("--additional_properties", metavar="PROPS", action="store",
                    type=str, default="", required=False,
                    help="Comma-separated list of additional properties to test (e.g., 'mass,type,v')")
parser.add_argument("--n_iter", metavar="N", action="store",
                    type=int, default=20, required=False,
                    help="Number of iterations for timing (default: 10)")

args = parser.parse_args()

# Parse additional properties
additional_props = []
if args.additional_properties:
    additional_props = [p.strip()
                        for p in args.additional_properties.split(',')]

# System setup
system = espressomd.System(box_l=[10.0, 10.0, 10.0])
system.time_step = 0.01
system.cell_system.skin = 0.5

n_proc = system.cell_system.get_state()["n_nodes"]
n_part = n_proc * args.particles_per_core

print(f"Benchmark Configuration:")
print(f"  Particles per core: {args.particles_per_core}")
print(f"  Total particles: {n_part}")
print(f"  MPI ranks: {n_proc}")
print(f"  Iterations: {args.n_iter}")
print(f"  Additional properties: {
      additional_props if additional_props else 'None'}")
print()

# Test 1: Bulk particle creation without providing IDs
print_section_header("Test 1: Bulk particle creation without IDs")


def create_particles_no_ids():
    system.part.clear()
    positions = np.random.random((n_part, 3)) * system.box_l
    system.part.add(pos=positions)


avg, std = benchmark_operation(create_particles_no_ids, args.n_iter)
print(f"Average: {avg:.6f} s ± {std:.6f} s ({n_part / avg:.0f} particles/s)")
print()

# Test 2: Bulk particle creation with linear IDs
print_section_header("Test 2: Bulk particle creation with linear IDs")


def create_particles_with_ids():
    system.part.clear()
    positions = np.random.random((n_part, 3)) * system.box_l
    ids = np.arange(n_part, dtype=int)
    system.part.add(pos=positions, id=ids)


avg, std = benchmark_operation(create_particles_with_ids, args.n_iter)
print(f"Average: {avg:.6f} s ± {std:.6f} s ({n_part / avg:.0f} particles/s)")
print()

# Setup particles for property access tests
system.part.clear()
positions = np.random.random((n_part, 3)) * system.box_l
system.part.add(pos=positions)

# Test 3: Single particle property access
print_section_header("Test 3: Single particle property access")

# Base properties to test
base_properties_read = ['pos', 'pos_folded', 'q', 'f']
base_properties_write = ['pos', 'q']

# Add additional properties and filter by availability
test_particle = system.part.by_id(0)
all_props_read = filter_available_properties(
    base_properties_read + additional_props, test_particle)
all_props_write = filter_available_properties(
    base_properties_write +
    [p for p in additional_props if p not in ['f', 'pos_folded']],
    test_particle)


def benchmark_property_access(properties, access_func, time_scale=1e6, time_unit="µs"):
    """Benchmark property access (read or write) and print results."""
    for prop in properties:
        avg, std = benchmark_operation(
            lambda: access_func(prop), args.n_iter)
        print(f"  {prop}: {avg * time_scale:.3f} {time_unit} ± {std *
              time_scale:.3f} {time_unit}, {1 / avg:.0f} particles/s")


# Test single particle read
print("\nSingle particle property READ:")
benchmark_property_access(
    all_props_read,
    lambda prop: getattr(test_particle, prop)
)

# Test single particle write
print("\nSingle particle property WRITE:")


def write_property(prop):
    test_value = getattr(test_particle, prop)
    setattr(test_particle, prop, test_value)


benchmark_property_access(all_props_write, write_property)
print()

# Test 4: Slice property access (re-instantiating slice)
print_section_header("Test 4: Slice property access (re-instantiating slice)")


def benchmark_slice_access(properties, access_func, description, time_scale=1e3, time_unit="ms"):
    """Benchmark slice property access and print results."""
    print(f"\n{description}:")
    for prop in properties:
        avg, std = benchmark_operation(
            lambda: access_func(prop), args.n_iter)
        print(f"  {prop}: {avg * time_scale:.3f} {time_unit} ± {std *
              time_scale:.3f} {time_unit} ({n_part / avg:.0f} particles/s)")


# Test slice read - re-instantiating slice on every access
benchmark_slice_access(
    all_props_read,
    lambda prop: getattr(system.part.all(), prop),
    "Slice property READ (re-instantiating slice)"
)

# Test slice write - re-instantiating slice on every access


def write_slice_property(prop):
    test_value = getattr(test_particle, prop)
    test_array = np.array([test_value] * n_part)
    setattr(system.part.all(), prop, test_array)


benchmark_slice_access(
    all_props_write,
    write_slice_property,
    "Slice property WRITE (re-instantiating slice)"
)
print()

# Test 5: Slice property access (cached slice)
print_section_header("Test 5: Slice property access (cached slice)")


def benchmark_cached_slice_access(properties, access_func, warmup_func, description):
    """Benchmark cached slice property access with warmup."""
    print(f"\n{description}:")
    for prop in properties:
        timings = []
        for _ in range(args.n_iter):
            particle_slice = system.part.all()
            # Warm-up access
            warmup_func(particle_slice, prop)
            # Timed access
            tick = time.time()
            access_func(particle_slice, prop)
            tock = time.time()
            timings.append(tock - tick)
        avg = np.mean(timings)
        std = np.std(timings)
        print(f"  {prop}: {avg * 1e3:.3f} ms ± {std *
              1e3:.3f} ms ({n_part / avg:.0f} particles/s)")


# Test slice read - cached slice
benchmark_cached_slice_access(
    all_props_read,
    lambda slice_obj, prop: getattr(slice_obj, prop),
    lambda slice_obj, prop: getattr(slice_obj, prop),
    "Slice property READ (cached slice)"
)

# Test slice write - cached slice


def write_cached_slice(slice_obj, prop):
    test_value = getattr(test_particle, prop)
    test_array = np.array([test_value] * n_part)
    setattr(slice_obj, prop, test_array)


benchmark_cached_slice_access(
    all_props_write,
    write_cached_slice,
    write_cached_slice,
    "Slice property WRITE (cached slice)"
)

print()
print_section_header("Benchmark completed")
