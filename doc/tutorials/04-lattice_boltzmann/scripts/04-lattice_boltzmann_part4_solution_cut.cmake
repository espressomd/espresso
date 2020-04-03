# This code splits 04-lattice_boltzmann_part4_solution.py such that it can be
# added to the tutorial without having two espressomd.System class instances
file(
  READ
  "${CMAKE_CURRENT_SOURCE_DIR}/scripts/04-lattice_boltzmann_part4_solution.py"
  script_full)
string(FIND "${script_full}" "# Extract fluid velocity along the x-axis"
            cut_position)
string(SUBSTRING "${script_full}" ${cut_position} -1 script_cut)
string(PREPEND script_cut "import matplotlib.pyplot as plt\n")
file(
  WRITE
  "${CMAKE_CURRENT_BINARY_DIR}/scripts/04-lattice_boltzmann_part4_solution_cut.py"
  "${script_cut}")
