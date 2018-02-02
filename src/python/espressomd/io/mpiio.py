from ..script_interface import PScriptInterface


class Mpiio(object):

    """MPI-IO object.

    Used to output particle data using MPI-IO to binary files.

    .. note::
        See the write and read methods for documentation.
    """

    def __init__(self):
        self._instance = PScriptInterface(
            "ScriptInterface::MPIIO::MPIIOScript")

    def write(self, prefix=None, positions=False, velocities=False, types=False, bonds=False):
        """MPI-IO write.

        Outputs binary data using MPI-IO to several files starting with prefix.
        Suffixes are:

            - head: Information about fields that are dumped,
            - pref: Information about processes: 1 int per process,
            - id: Particle ids: 1 int per particle,
            - pos: Position information (if dumped): 3 doubles per particle,
            - vel: Velocity information (if dumped): 3 doubles per particle,
            - typ: Type information (if dumped): 1 int per particle,
            - bond: Bond information (if dumped): variable amount of data,
            - boff: Bond offset information (if bonds are dumped): 1 int per particle.

        .. note::
            Do not read the files on a machine with a different architecture!

        Parameters
        ----------
        prefix : :obj:`str`
            Common prefix for the filenames.
        positions : :obj:`bool`, optional
            Indicates if positions should be dumped.
        velocities : :obj:`bool`, optional
            Indicates if velocities should be dumped.
        types : :obj:`bool`, optional
            Indicates if types should be dumped.
        bonds : :obj:`bool`, optional
            Indicates if bonds should be dumped.

        Raises
        ------
        ValueError
            If no prefix was given or none of the output fields are chosen.
        """

        if prefix is None:
            raise ValueError(
                "Need to supply output prefix via 'prefix' kwarg.")
        if not positions and not velocities and not types and not bonds:
            raise ValueError("No output fields chosen.")

        self._instance.call_method(
            "write", prefix=prefix, pos=positions, vel=velocities, typ=types, bond=bonds)

    def read(self, prefix=None, positions=False, velocities=False, types=False, bonds=False):
        """MPI-IO read.

        This function reads data dumped by `write`. See the write documentation
        for details.

        .. note::
            The files must be read on the same number of processes that wrote
            the data. The data must be read on a machine with the same
            architecture (otherwise, this might silently fail).
        """
        if prefix is None:
            raise ValueError(
                "Need to supply output prefix via 'prefix' kwarg.")
        if not positions and not velocities and not types and not bonds:
            raise ValueError("No output fields chosen.")

        self._instance.call_method(
            "read", prefix=prefix, pos=positions, vel=velocities, typ=types, bond=bonds)

mpiio = Mpiio()
