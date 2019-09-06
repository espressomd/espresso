# Copyright (C) 2010-2018 The ESPResSo project
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
from .script_interface import ScriptInterfaceHelper, script_interface_register
from .particle_data import ParticleHandle, ParticleSlice


@script_interface_register
class Cluster(ScriptInterfaceHelper):

    """Class representing a cluster of particles

    Methods
    -------
    particle_ids():
        Returns list of particle ids in the cluster

    particles():
        Returns an instance of ParticleSlice containing the particles in the cluster

    size():
        Returns the number of particles in the cluster

    center_of_mass():
        center of mass of the cluster

    longest_distance():
        Longest distance between any combination of two particles in the cluster

    fractal_dimension(dr=None):
        estimates the cluster's fractal dimension by fitting the number of particles
        :math:`n` in spheres of growing radius around the center of mass
        to :math:`c*r_g^d`, where :math:`r_g` is the radius of gyration of the particles
        within the sphere, and :math:`d` is the fractal dimension.
        `dr`: Minimum increment for the radius of the spheres.
        Return value: (fractal_dimension, mean_square_residual)
    """
    _so_name = "ClusterAnalysis::Cluster"
    _so_bind_methods = ("particle_ids", "size", "longest_distance",
                        "radius_of_gyration", "fractal_dimension", "center_of_mass")

    _so_creation_policy = "LOCAL"

    def particles(self):
        return ParticleSlice(self.particle_ids())


@script_interface_register
class ClusterStructure(ScriptInterfaceHelper):

    """Cluster structure of a simulation system, and access to cluster analysis

    Attributes
    ----------
    pair_criterion: classes derived from ``_PairCriterion``
        Criterion to decide whether two particles are neighbors.

    clusters: behaves like a read-only dictionary
        Access to individual clusters in the cluster structure either via
        cluster[i], where i is a (non-consecutive) integer cluster id
        or via iteration::

            for pair in clusters:

        where pair contains the numeric id and the corresponding cluster object.

    """
    _so_name = "ClusterAnalysis::ClusterStructure"
    _so_creation_policy = "LOCAL"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._clusters = Clusters(self)

    def run_for_all_pairs(self):
        """
        Runs the cluster analysis, considering all pairs of particles in the system

        """
        return self.call_method("run_for_all_pairs")

    def run_for_bonded_particles(self):
        """
        Runs the cluster analysis, considering only pairs of particles connected by a pair-bond.

        """
        return self.call_method("run_for_bonded_particles")

    def clear(self):
        """
        Clears the cluster structure.

        """
        return self.call_method("clear")

    def cluster_ids(self):
        """
        Returns a list of all cluster ids of the clusters in the structure.

        """
        return self.call_method("cluster_ids")

    def cid_for_particle(self, p):
        """Returns cluster id for the particle (passed as ParticleHandle or particle id)"""
        if isinstance(p, ParticleHandle):
            return self.call_method("cid_for_particle", pid=p.id)
        if isinstance(p, int):
            return self.call_method("cid_for_particle", pid=p)
        else:
            raise TypeError(
                "The particle has to be passed as instance of Particle handle or as an integer particle id")

    @property
    def clusters(self):
        """Gives access to the clusters in the cluster structure via an instance of :any:`Clusters`."""
        return self._clusters


class Clusters:

    """Access to the clusters in the cluster structure.

       Access is as follows:

       * Number of clusters: len(clusters)
       * Access a cluster via its id: clusters[id]
       * Iterate over clusters::

            for c in clusters:

         where c will be a tuple containing the cluster id and the cluster object
     """

    def __init__(self, cluster_structure):
        self.cluster_structure = cluster_structure

    def __getitem__(self, cluster_id):
        return self.cluster_structure.call_method("get_cluster", id=cluster_id)

    def __iter__(self):
        for cid in self.cluster_structure.cluster_ids():
            yield (cid, self.cluster_structure.call_method("get_cluster", id=cid))

    def __len__(self):
        return self.cluster_structure.call_method("n_clusters")
