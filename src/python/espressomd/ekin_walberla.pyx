include "myconfig.pxi"

import itertools
import functools

from .actors cimport Actor
from .utils import array_locked, to_char_pointer, check_type_or_throw_except
from .utils cimport Vector3i, Vector3d, make_array_locked, create_nparray_from_double_array

from .lb import _construct

import numpy as np
import os

# TODO: figure out duplication between LB and EKin slices
# TODO: consistent "ekin"/"ek" prefix
# TODO: figure out how to remove duplication in the MPI-calls between ekin/lb


# TODO: boundaries?
# TODO: vtk writer

IF EK_WALBERLA:

    # TODO: this is a full duplicate of LB
    import lxml.etree

    class VTKRegistry:

        def __init__(self):
            self.map = {}
            self.collisions = {}

        def _register_vtk_object(self, vtk_uid, vtk_obj):
            self.map[vtk_uid] = vtk_obj
            self.collisions[os.path.abspath(vtk_uid)] = vtk_uid

        def __getstate__(self):
            return self.map

        def __setstate__(self, active_vtk_objects):
            self.map = {}
            self.collisions = {}
            for vtk_uid, vtk_obj in active_vtk_objects.items():
                self.map[vtk_uid] = vtk_obj
                self.collisions[os.path.abspath(vtk_uid)] = vtk_uid

    _vtk_registry = VTKRegistry()

    class VTKOutput:
        """
        VTK callback.
        """
        observable2enum = {'density': < int > output_vtk_density}

        def __init__(self, vtk_uid, identifier, observables, delta_N,
                     base_folder, prefix):
            observables = set(observables)
            flag = sum(self.observable2enum[obs] for obs in observables)
            self._params = {
                'vtk_uid': vtk_uid, 'identifier': identifier, 'prefix': prefix,
                'delta_N': delta_N, 'base_folder': base_folder, 'flag': flag,
                'enabled': True, 'initial_count': 0, 'observables': observables}
            self._set_params_in_es_core()

        def __getstate__(self):
            odict = self._params.copy()
            # get initial execution counter
            pvd = os.path.abspath(self._params['vtk_uid']) + '.pvd'
            if os.path.isfile(pvd):
                tree = lxml.etree.parse(pvd)
                nodes = tree.xpath('/VTKFile/Collection/DataSet')
                if nodes:
                    odict['initial_count'] = int(
                        nodes[-1].attrib['timestep']) + 1
            return odict

        def __setstate__(self, params):
            self._params = params
            self._set_params_in_es_core()

        def _set_params_in_es_core(self):
            _vtk_registry._register_vtk_object(self._params['vtk_uid'], self)
            create_vtk(self._params['delta_N'],
                       self._params['initial_count'],
                       self._params['flag'],
                       to_char_pointer(self._params['identifier']),
                       to_char_pointer(self._params['base_folder']),
                       to_char_pointer(self._params['prefix']))
            if not self._params['enabled']:
                self.disable()

        def write(self):
            raise NotImplementedError()

        def disable(self):
            raise NotImplementedError()

        def enable(self):
            raise NotImplementedError()

    class VTKOutputAutomatic(VTKOutput):
        """
        Automatic VTK callback. Can be disabled at any time and re-enabled later.

        Please note that the internal VTK counter is no longer incremented
        when the callback is disabled, which means the number of LB steps
        between two frames will not always be an integer multiple of delta_N.
        """

        def __repr__(self):
            return "<{}.{}: writes to '{}' every {} LB step{}{}>".format(
                self.__class__.__module__, self.__class__.__name__,
                self._params['vtk_uid'], self._params['delta_N'],
                '' if self._params['delta_N'] == 1 else 's',
                '' if self._params['enabled'] else ' (disabled)')

        def write(self):
            raise RuntimeError(
                'Automatic VTK callbacks writes automaticall, cannot be triggered manually')

        def enable(self):
            switch_vtk(to_char_pointer(self._params['vtk_uid']), 1)
            self._params['enabled'] = True

        def disable(self):
            switch_vtk(to_char_pointer(self._params['vtk_uid']), 0)
            self._params['enabled'] = False

    class VTKOutputManual(VTKOutput):
        """
        Manual VTK callback. Can be called at any time to take a snapshot
        of the current state of the LB fluid.
        """

        def __repr__(self):
            return "<{}.{}: writes to '{}' on demand>".format(
                self.__class__.__module__, self.__class__.__name__,
                self._params['vtk_uid'])

        def write(self):
            write_vtk(to_char_pointer(self._params['vtk_uid']))

        def disable(self):
            raise RuntimeError('Manual VTK callbacks cannot be disabled')

        def enable(self):
            raise RuntimeError('Manual VTK callbacks cannot be enabled')

    cdef class EKinWalberla(Actor):
        def validate_params(self):
            pass

        def valid_keys(self):
            return {"diffusion", "kT", "dens", "tau"}

        def required_keys(self):
            return {"diffusion", "kT", "dens", "tau"}

        def default_params(self):
            return {}

        def _set_params_in_es_core(self):
            pass

        def _get_params_from_es_core(self):
            self._params['diffusion'] = self.diffusion
            self._params["kT"] = self.kT
            self._params["tau"] = self.tau

            return self._params

        def __reduce__(self):
            return _construct, (self.__class__, self._params), None

        def __getitem__(self, key):
            cdef Vector3i shape
            if isinstance(key, (tuple, list, np.ndarray)):
                if len(key) == 3:
                    if any(isinstance(typ, slice) for typ in key):
                        shape = get_shape()
                        return EKinSlice(key, (shape[0], shape[1], shape[2]))
                    else:
                        return EKinRoutines(np.array(key))
            else:
                raise Exception(
                    f"{key} is not a valid key. Should be a point on the nodegrid e.g. ek[0,0,0], or a slice")

        def save_checkpoint(self, path, binary):
            # tmp_path = path + ".__tmp__"
            # lb_lbfluid_save_checkpoint(utils.to_char_pointer(tmp_path), binary)
            # os.rename(tmp_path, path)
            pass

        def load_checkpoint(self, path, binary):
            # lb_lbfluid_load_checkpoint(utils.to_char_pointer(path), binary)
            pass

        def _activate_method(self):
            mpi_init_ekin_walberla(
                self._params["diffusion"],
                self._params["kT"],
                self._params["dens"],
                self._params["tau"])

        def _deactivate_method(self):
            mpi_destruct_ekin_walberla()

        property kT:
            def __get__(self):
                return get_kT()

            def __set__(self, kT):
                set_kT(kT)

        property diffusion:
            def __get__(self):
                return get_diffusion()

            def __set__(self, diffusion):
                set_diffusion(diffusion)

        property tau:
            def __get__(self):
                return get_tau()

        def nodes(self):
            """Provides a generator for iterating over all lb nodes"""

            shape = self.shape
            for i, j, k in itertools.product(
                    range(shape[0]), range(shape[1]), range(shape[2])):
                yield self[i, j, k]

        def add_vtk_writer(self, identifier, observables, delta_N=0,
                           base_folder='vtk_out', prefix='simulation_step'):
            """
            Create a VTK observable.

            Files are written to ``<base_folder>/<identifier>/<prefix>_*.vtu``.
            Summary is written to ``<base_folder>/<identifier>.pvd``.

            Parameters
            ----------
            identifier : :obj:`str`
                Name of the VTK dataset.
            observables : :obj:`list`, \{'density', 'velocity_vector', 'pressure_tensor'\}
                List of observables to write to the VTK files.
            delta_N : :obj:`int`
                Write frequency, if 0 write a single frame (default),
                otherwise add a callback to write every ``delta_N`` LB steps
                to a new file.
            base_folder : :obj:`str`
                Path to the output VTK folder.
            prefix : :obj:`str`
                Prefix for VTK files.
            """
            # sanity checks
            check_type_or_throw_except(
                delta_N, 1, int, "delta_N must be 1 integer")
            assert delta_N >= 0, 'delta_N must be a positive integer'
            check_type_or_throw_except(
                identifier, 1, str, "identifier must be 1 string")
            assert os.path.sep not in identifier, 'identifier must be a string, not a filepath'
            if isinstance(observables, str):
                observables = [observables]
            for obs in observables:
                if obs not in VTKOutput.observable2enum:
                    raise ValueError('Unknown VTK observable ' + obs)
            vtk_uid = base_folder + '/' + identifier
            vtk_path = os.path.abspath(vtk_uid)
            if vtk_path in _vtk_registry.collisions:
                raise RuntimeError(
                    'VTK identifier "{}" would overwrite files written by VTK identifier "{}"'.format(
                        vtk_uid, _vtk_registry.collisions[vtk_path]))
            args = (
                vtk_uid,
                identifier,
                observables,
                delta_N,
                base_folder,
                prefix)
            if delta_N:
                obj = VTKOutputAutomatic(*args)
            else:
                obj = VTKOutputManual(*args)
            return obj


class EKinSlice:

    def __init__(self, key, shape):
        self.x_indices, self.y_indices, self.z_indices = self.get_indices(
            key, shape[0], shape[1], shape[2])

    def get_indices(self, key, shape_x, shape_y, shape_z):
        x_indices = np.atleast_1d(np.arange(shape_x)[key[0]])
        y_indices = np.atleast_1d(np.arange(shape_y)[key[1]])
        z_indices = np.atleast_1d(np.arange(shape_z)[key[2]])
        return x_indices, y_indices, z_indices

    def get_values(self, x_indices, y_indices, z_indices, prop_name):
        shape_res = np.shape(
            getattr(EKinRoutines(np.array([0, 0, 0])), prop_name))
        res = np.zeros(
            (x_indices.size,
             y_indices.size,
             z_indices.size,
             *shape_res))
        for i, x in enumerate(x_indices):
            for j, y in enumerate(y_indices):
                for k, z in enumerate(z_indices):
                    res[i, j, k] = getattr(EKinRoutines(
                        np.array([x, y, z])), prop_name)
        if shape_res == (1,):
            res = np.squeeze(res, axis=-1)
        return array_locked(res)

    def set_values(self, x_indices, y_indices, z_indices, prop_name, value):
        for i, x in enumerate(x_indices):
            for j, y in enumerate(y_indices):
                for k, z in enumerate(z_indices):
                    setattr(EKinRoutines(
                        np.array([x, y, z])), prop_name, value[i, j, k])


cdef class EKinRoutines:
    def __init__(self, key):
        check_type_or_throw_except(
            key, 3, int, "The index of an lb fluid node consists of three integers.")
        self.node[0] = key[0]
        self.node[1] = key[1]
        self.node[2] = key[2]
        if not node_is_index_valid(self.node):
            raise ValueError("EKin node index out of bounds")

    property index:
        def __get__(self):
            return (self.node[0], self.node[1], self.node[2])

    property density:
        def __get__(self):
            return get_density(self.node)

        def __set__(self, value):
            set_density(self.node, value)

    property is_boundary:
        def __get__(self):
            return get_is_boundary(self.node)


def _add_ek_slice_properties():
    """
    Automatically add all of EKinRoutines's properties to LBSlice.

    """

    def set_attribute(ek_slice, value, attribute):
        """
        Setter function that sets attribute on every member of lb_slice.
        If values contains only one element, all members are set to it.

        """

        indices = [ek_slice.x_indices, ek_slice.y_indices, ek_slice.z_indices]
        N = [len(x) for x in indices]

        if N[0] * N[1] * N[2] == 0:
            raise AttributeError("Cannot set properties of an empty LBSlice")

        value = np.copy(value)
        attribute_shape = ek_slice.get_values(
            *np.zeros((3, 1), dtype=int), attribute).shape[3:]
        target_shape = (*N, *attribute_shape)

        # broadcast if only one element was provided
        if value.shape == attribute_shape:
            value = np.ones(target_shape) * value

        if value.shape != target_shape:
            raise ValueError(
                f"Input-dimensions of {attribute} array {value.shape} does not match slice dimensions {target_shape}.")

        ek_slice.set_values(*indices, attribute, value)

    def get_attribute(ek_slice, attribute):
        """
        Getter function that copies attribute from every member of
        lb_slice into an array (if possible).

        """

        indices = [ek_slice.x_indices, ek_slice.y_indices, ek_slice.z_indices]
        N = [len(x) for x in indices]

        if N[0] * N[1] * N[2] == 0:
            return np.empty(0, dtype=type(None))

        return ek_slice.get_values(*indices, attribute)

    for attribute_name in dir(EKinRoutines):
        if attribute_name in dir(EKinSlice) or not isinstance(
                getattr(EKinRoutines, attribute_name), type(EKinRoutines.density)):
            continue

        # synthesize a new property
        new_property = property(
            functools.partial(get_attribute, attribute=attribute_name),
            functools.partial(set_attribute, attribute=attribute_name),
            doc=getattr(EKinRoutines, attribute_name).__doc__ or f'{attribute_name} for a slice')
        # attach the property to LBSlice
        setattr(EKinSlice, attribute_name, new_property)


_add_ek_slice_properties()
