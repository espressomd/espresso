# Guidelines for Coding Agents Working on ESPResSo

## Purpose of the Software

ESPResSo is a framework for simulations in soft matter science, statistical physics, process engineering, and related fields.

**Key Algorithms**:
- Molecular Dynamics (MD)
- Monte Carlo methods
- P3M electrostatics (Particle-Particle-Particle-Mesh)
- Lattice-Boltzmann hydrodynamics
- Electrokinetics via diffusion-advection-reaction solver coupled to Poisson solver

## Usage

Simulations are controlled through Python scripts via the `espressomd` module.

**Usage Examples**:
- `samples/` - Sample simulations demonstrating various features
- `testsuite/python/` - Comprehensive test suite with usage patterns
- `doc/tutorials/` - Step-by-step tutorials for learning ESPResSo

## Repository Structure

```
/src/core/                  # Simulation core written in C++
/src/python/                # Python bindings for espressomd module
/src/script_interface/      # Interface layer between Python and MPI-parallel C++ objects
/src/walberla_bridge/       # Interface library for lattice-based methods (LB, EK)
/src/utils/                 # Utility functions and helper classes
/testsuite/python/          # Python-based test suite
/doc/                       # Documentation source files
/samples/                   # Example simulation scripts
/build/                     # Build directory (created by user)
```

## Building

**Build System**: CMake-based configuration and compilation

**Standard Build Process**:
```bash
cd build
make -j$(nproc)
```

**Important Notes**:
- Builds should be done in the `build/` directory
- **Do not run cmake** unless:
  - The build folder does not yet contain CMake build files, OR
  - The user explicitly asks you to do so
- Use `$(nproc)` for optimal parallel compilation

## Testing

### Python Tests

**Basic Testing**:
```bash
cd build
./pypresso ../testsuite/python/test_script.py
```

**Parallel Testing** (MPI + OpenMP):
```bash
OMP_NUM_THREADS=4 mpirun -n 4 ./pypresso ../testsuite/python/test_script.py
```

### C++ Unit Tests

**Location**:
- `src/core/unit_tests/` - Core functionality unit tests
- `src/walberla_bridge/tests/` - Lattice-Boltzmann and EK tests

**Running C++ Tests**:
1. Build the test executable via CMake
2. Run directly from the build directory:
   ```bash
   cd build
   ./src/core/unit_tests/test_name
   ```

## Code Formatting

**Python Files**:
```bash
maintainer/format/autopep8.sh -i file.py
```

**C++ Files**:
```bash
maintainer/format/clang-format.sh -i file.cpp
```

**Important**: If formatters complain about missing Python packages or wrong versions, **stop and ask the user to load the environment**.

## Important Header Files in ESPResSo

### Core Simulation Framework

**System Management**:
- `src/core/system/System.hpp` - Main system class managing all simulation components
- `src/core/BoxGeometry.hpp` - Simulation box geometry and periodic boundary conditions

**Particle Data & Geometry**:
- `src/core/Particle.hpp` - Core particle data structure (position, velocity, force, properties)
- `src/core/ParticleRange.hpp` - Iterator interface for particle collections
- `src/core/cell_system/CellStructure.hpp` - Cell system interface for particle storage and spatial decomposition of the system

**Integration & Forces**:
- `src/core/integrate.hpp` - Main integration loop
- `src/core/forces.hpp` - Force calculation framework

### Interactions

**Non-bonded Interactions**:
- `src/core/nonbonded_interactions/nonbonded_interaction_data.hpp` - Parameter storage

**Bonded Interactions**:
- `src/core/bonded_interactions/bonded_interaction_data.hpp` - Parameter storage

**Electrostatics**:
- `src/core/electrostatics/coulomb.hpp` - Coulomb interaction interface
- `src/core/electrostatics/p3m.hpp` - P3M algorithm for long-range Coulomb interactions

### Thermostats & Observables

**Thermostats**:
- `src/core/thermostat.hpp` - Main thermostat interface (Langevin, DPD, Brownian, NPT)
- `src/core/thermostats/` - Thermostat implementations

**Observables & Statistics**:
- `src/core/Observable_stat.hpp` - Observable statistics collection

### Utilities & Communication

**Communication**:
- `src/core/communication.hpp` - MPI communication management

**Mathematical Utilities**:
- `src/utils/include/utils/Vector.hpp` - 3D vector math operations

## Further Information

**Documentation**:
- User Guide: `doc/sphinx/` (Sphinx-based documentation)
- Developer Wiki: https://github.com/espressomd/espresso/wiki
- Tutorials: `doc/tutorials/`

**Getting Help**:
- Check the developer wiki for architectural overviews
- Consult existing tests in `testsuite/python/` for usage patterns
- Review sample scripts in `samples/` for common simulation setups

## Script Interface

The script interface handles the representation of Python objects as C++ objects in the parallel simulation core. It provides the bridge between Python-level user code and MPI-parallel C++ simulation objects.

### Architecture Overview

**Purpose**: Enable Python control of parallel C++ simulation objects with automatic parameter synchronization across MPI ranks.

**Key Components**:
1. **Python Layer**: Cython-based wrappers that expose C++ objects to Python
2. **Interface Layer**: C++ classes that handle object lifecycle and parameter management
3. **Core Layer**: Actual simulation objects in the parallel core

### Entry Points and Important Files

**Python Side**:
- `src/python/espressomd/script_interface.pyx` - Core infrastructure
  - `ScriptInterfaceHelper` class - Base class for Python-side script interface objects
  - `script_interface_register` decorator - Registers C++ classes for Python access

**C++ Side**:
- `src/script_interface/ObjectHandle.hpp` - Base class for all C++ script interface objects
- `src/script_interface/auto_parameters/AutoParameters.hpp` - Automatic parameter handling
  - Provides `AutoParameters` template for objects with parameters
  - Used for most script interface classes to reduce boilerplate

### Implementation Patterns

#### Example 1: Simple Objects (Shapes)

**Use Case**: Individual objects with parameters (e.g., geometric shapes)

**Python Implementation**:
- `src/python/espressomd/shapes.py` - Python shape classes

**C++ Implementation**:
- `src/script_interface/shapes/` - C++ shape interface classes
  - `Cylinder.hpp`, `Ellipsoid.hpp`, etc. - Individual shape implementations
  - Each inherits from `AutoParameters<Shape, ...>` for automatic parameter handling

**Pattern**: One Python class maps to one C++ interface class, which wraps a core simulation object.

#### Example 2: Container Objects (Constraints)

**Use Case**: Collections that manage multiple objects (e.g., constraint container)

**Python Implementation**:
- `src/python/espressomd/constraints.py` - Constraint container and constraint classes

**C++ Implementation**:
- `src/script_interface/constraints/Constraints.hpp` - Container interface class
  - Manages collection of constraint objects
  - Provides add/remove functionality
  - Synchronizes across MPI ranks

**Pattern**: Container classes manage collections of script interface objects with automatic MPI synchronization.

### Best Practices

**When Creating New Script Interface Objects**:
1. Inherit from `AutoParameters` for automatic parameter handling
2. Register parameters in the constructor using `add_parameters()`
3. Use `script_interface_register` decorator on the Python side
4. Follow existing patterns in similar objects (shapes, constraints, etc.)
5. Ensure proper MPI synchronization for parallel operations

**Common Parameter Types**:
- Scalar values (int, double, bool)
- Vectors (`Utils::Vector3d`, `std::vector<T>`)
- References to other script interface objects
- Variant types for flexible parameters

