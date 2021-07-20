include "myconfig.pxi"

import itertools
import functools

from .actors cimport Actor
from .utils import array_locked
from .utils cimport Vector3i, Vector3d, make_array_locked, create_nparray_from_double_array
from .utils import check_type_or_throw_except

from .lb import _construct

import numpy as np

# TODO: figure out duplication between LB and EKin slices
# TODO: consistent "ekin"/"ek" prefix
# TODO: figure out how to remove duplication in the MPI-calls between ekin/lb


# TODO: boundaries?
# TODO: vtk writer

IF EK_WALBERLA:
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
                    "%s is not a valid key. Should be a point on the nodegrid e.g. lbf[0,0,0], or a slice" % key)

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
