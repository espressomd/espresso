# Copyright (C) 2010-2019 The ESPResSo project
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
from . import interactions
from .__init__ import has_features


class DrudeHelpers:
    """
    Provides helper functions to assist in the creation of Drude particles.

    """

    def __init__(self):
        # Dict with Drude type infos
        self.drude_dict = {}
        # Lists with unique Drude and core types
        self.core_type_list = []
        self.drude_type_list = []
        # Get core id from Drude id
        self.core_id_from_drude_id = {}
        # Drude IDs
        self.drude_id_list = []

    def add_drude_particle_to_core(self, system, harmonic_bond,
                                   thermalized_bond, p_core, type_drude, alpha,
                                   mass_drude, coulomb_prefactor,
                                   thole_damping=2.6, verbose=False):
        """
        Adds a Drude particle with specified type and mass to the system and
        returns its `espressomd.particle_data.ParticleHandle`.
        Checks if different Drude particles have different types.
        Collects types/charges/polarizations/Thole factors for intramolecular
        core-Drude short-range exclusion and Thole interaction.

        Parameters
        ----------

        system : :class:`espressomd.system.System`
        harmonic_bond: :class:`espressomd.interactions.HarmonicBond`
            Add this harmonic bond to between Drude particle and core
        thermalized_bond: :class:`espressomd.interactions.ThermalizedBond`
            Add this thermalized_bond to between Drude particle and core
        p_core: :class:`espressomd.particle_data.ParticleHandle`
            The existing core particle
        type_drude: :obj:`int`
            The type of the newly created Drude particle
        alpha : :obj:`float`
            The polarizability in units of inverse volume. Related to the charge
            of the Drude particle.
        mass_drude : :obj:`float`
            The mass of the newly created Drude particle
        coulomb_prefactor : :obj:`float`
            Required to calculate the charge of the Drude particle.
        thole_damping : :obj:`float`
            Thole damping factor of the Drude pair. Comes to effect if
            :meth:`add_all_thole()` method is used.
        verbose : :obj:`bool`
            Turns on verbosity.

        Returns
        -------
        :class:`espressomd.particle_data.ParticleHandle`
            The created Drude Particle.

        """

        k = harmonic_bond.params["k"]
        q_drude = -1.0 * pow(k * alpha / coulomb_prefactor, 0.5)

        if has_features("PARTICLE_ANISOTROPY"):
            gamma_off = [0.0, 0.0, 0.0]
        else:
            gamma_off = 0.0

        drude_part = system.part.add(pos=p_core.pos, type=type_drude,
                                     q=q_drude, mass=mass_drude, gamma=gamma_off)
        id_drude = drude_part.id

        if verbose:
            print(
                f"Adding to core {p_core.id}   pol {alpha}   core charge {p_core.q} -> {p_core.q - q_drude}   drude charge {q_drude}")

        p_core.q -= q_drude
        p_core.mass -= mass_drude
        p_core.add_bond((harmonic_bond, id_drude))
        p_core.add_bond((thermalized_bond, id_drude))

        p_core.gamma = gamma_off

        if type_drude in self.drude_dict and not (
                self.drude_dict[type_drude]["q"] == q_drude and
                self.drude_dict[type_drude]["thole_damping"] == thole_damping):
            raise Exception(
                "Drude particles with different drude charges have to have different types for THOLE")

        self.core_id_from_drude_id[id_drude] = p_core.id

        # Add new Thole non-bonded interaction for D-D, D-C, C-C for all existing
        # Drude types if this type is seen for the first time
        if type_drude not in self.drude_dict:

            # Bookkeeping of q, alphas and damping parameter
            self.drude_dict[type_drude] = {}
            self.drude_dict[type_drude]["q"] = q_drude
            self.drude_dict[type_drude]["qc"] = p_core.q
            self.drude_dict[type_drude]["alpha"] = alpha
            self.drude_dict[type_drude]["thole_damping"] = thole_damping
            self.drude_dict[type_drude]["core_type"] = p_core.type
            # Save same information to get access to the parameters via core
            # types
            self.drude_dict[p_core.type] = {}
            self.drude_dict[p_core.type]["q"] = -q_drude
            self.drude_dict[p_core.type]["qc"] = p_core.q
            self.drude_dict[p_core.type]["alpha"] = alpha
            self.drude_dict[p_core.type]["thole_damping"] = thole_damping
            self.drude_dict[p_core.type]["drude_type"] = type_drude

        # Collect unique Drude types
        if type_drude not in self.drude_type_list:
            self.drude_type_list.append(type_drude)

        # Collect unique core types
        if p_core.type not in self.core_type_list:
            self.core_type_list.append(p_core.type)

        # Collect unique Drude ids
        if id_drude not in self.drude_id_list:
            self.drude_id_list.append(id_drude)

        return drude_part

    def add_thole_pair_damping(self, system, t1, t2, verbose=False):
        """
        Calculates mixed Thole factors depending on Thole damping and polarization.
        Adds non-bonded Thole interactions to the system.

        Parameters
        ----------

        system : :class:`espressomd.system.System`
        t1 : :obj:`int`
            Type 1
        t2 : :obj:`int`
            Type 2
        verbose : :obj:`bool`
            Turns on verbosity.

        """

        qq = self.drude_dict[t1]["q"] * self.drude_dict[t2]["q"]
        s = 0.5 * (self.drude_dict[t1]["thole_damping"] + self.drude_dict[t2]["thole_damping"]
                   ) / (self.drude_dict[t1]["alpha"] * self.drude_dict[t2]["alpha"])**(1.0 / 6.0)

        system.non_bonded_inter[t1, t2].thole.set_params(
            scaling_coeff=s, q1q2=qq)

        if verbose:
            print(
                f"Added THOLE non-bonded interaction for types {t1} <-> {t2} S {s} q1q2 {qq}")

    def add_all_thole(self, system, verbose=False):
        """
        Calls :meth:`add_thole_pair_damping()` for all necessary combinations to
        create the interactions.

        Parameters
        ----------

        system : :class:`espressomd.system.System`
        verbose : :obj:`bool`
            Turns on verbosity.

        """

        # Drude <-> Drude
        for i in range(len(self.drude_type_list)):
            for j in range(i, len(self.drude_type_list)):
                self.add_thole_pair_damping(
                    system, self.drude_type_list[i], self.drude_type_list[j], verbose)
        # core <-> core
        for i in range(len(self.core_type_list)):
            for j in range(i, len(self.core_type_list)):
                self.add_thole_pair_damping(
                    system, self.core_type_list[i], self.core_type_list[j], verbose)
        # Drude <-> core
        for i in self.drude_type_list:
            for j in self.core_type_list:
                self.add_thole_pair_damping(system, i, j, verbose)

    def setup_and_add_drude_exclusion_bonds(self, system, verbose=False):
        """
        Creates electrostatic short-range exclusion bonds for global exclusion
        between Drude particles and core charges and adds the bonds to the cores.
        Has to be called once after all Drude particles have been created.

        Parameters
        ----------

        system : :class:`espressomd.system.System`
        verbose: :obj:`bool`
            Turns on verbosity.

        """

        # All Drude types need...
        for td in self.drude_type_list:

            # ...exclusions with core
            qd = self.drude_dict[td]["q"]  # Drude charge
            qc = self.drude_dict[td]["qc"]  # Core charge
            subtr_sr_bond = interactions.BondedCoulombSRBond(
                q1q2=-qd * qc)
            system.bonded_inter.add(subtr_sr_bond)
            self.drude_dict[td]["subtr_sr_bonds_drude-core"] = subtr_sr_bond
            if verbose:
                print(
                    f"Added drude-core SR exclusion bond {subtr_sr_bond} for drude {qd} <-> core {qc} to system")

        for drude_id in self.drude_id_list:
            core_id = self.core_id_from_drude_id[drude_id]
            pd = system.part.by_id(drude_id)
            pc = system.part.by_id(core_id)
            bond = self.drude_dict[pd.type]["subtr_sr_bonds_drude-core"]
            pc.add_bond((bond, drude_id))
            if verbose:
                print(
                    f"Added drude-core SR bond {bond} between ids {drude_id} and {core_id}")

    def setup_intramol_exclusion_bonds(self, system, mol_drude_types, mol_core_types,
                                       mol_core_partial_charges, verbose=False):
        """
        Creates electrostatic short-range exclusion bonds for intramolecular exclusion
        between Drude particles and partial charges of the cores. Has to be called once
        after all Drude particles have been created.

        Parameters
        ----------

        system : :class:`espressomd.system.System`
        mol_drude_types :
            List of types of Drude particles within the molecule
        mol_core_types :
            List of types of core particles within the molecule
        mol_core_partial_charges :
            List of partial charges of core particles within the molecule
        verbose : :obj:`bool`
            Turns on verbosity.

        """

        # All Drude types need...
        for td in mol_drude_types:
            self.drude_dict[td]["subtr_sr_bonds_intramol"] = {}

            # ... sr exclusion bond with other partial core charges...
            for tc, qp in zip(mol_core_types, mol_core_partial_charges):
                # ...excluding the Drude core partner
                if self.drude_dict[td]["core_type"] != tc:
                    qd = self.drude_dict[td]["q"]  # Drude charge
                    subtr_sr_bond = interactions.BondedCoulombSRBond(
                        q1q2=-qd * qp)
                    system.bonded_inter.add(subtr_sr_bond)
                    self.drude_dict[td]["subtr_sr_bonds_intramol"][
                        tc] = subtr_sr_bond
                    if verbose:
                        print(
                            f"Added intramolecular exclusion {subtr_sr_bond} for drude {qd} <-> core {qp} to system")

    def add_intramol_exclusion_bonds(
            self, drude_parts, core_parts, verbose=False):
        """
        Applies electrostatic short-range exclusion bonds for the given ids.
        Has to be applied for all molecules.

        Parameters
        ----------

        system : :class:`espressomd.system.System`
        drude_parts :
            List of Drude particles within a molecule.
        core_parts :
            List of core particles within a molecule.
        verbose : :obj:`bool`
            Turns on verbosity.

        """

        for drude_part in drude_parts:
            drude_id = drude_part.id
            for core_part in core_parts:
                core_id = core_part.id
                if self.core_id_from_drude_id[drude_id] != core_id:
                    pd = drude_part
                    pc = core_part
                    bond = self.drude_dict[pd.type][
                        "subtr_sr_bonds_intramol"][pc.type]
                    pd.add_bond((bond, core_id))
                    if verbose:
                        print(
                            f"Added subtr_sr bond {bond} between ids {drude_id} and {core_id}")
