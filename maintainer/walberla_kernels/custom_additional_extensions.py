#
# Copyright (C) 2022-2023 The ESPResSo project
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

import pathlib

import jinja2
import numpy as np
import pystencils as ps
import pystencils_walberla
import sympy as sp


class Dirichlet_Custom(ps.boundaries.Dirichlet):
    inner_or_boundary = False
    single_link = False  # this is the hacky solution

    def __init__(self, value, name=None, data_type="double"):
        super().__init__(value=value, name=name)
        self.data_type = data_type

    @property
    def additional_data(self):
        if callable(self._value):
            return [('value', ps.typing.BasicType(self.data_type))]
        else:
            return []


class Flux(ps.boundaries.boundaryconditions.Boundary):
    inner_or_boundary = True  # call the boundary condition with the fluid cell
    single_link = False  # needs to be called for all directional fluxes

    def __init__(self, stencil, value=None, dim=None, data_type='double'):
        self.stencil = stencil
        self.value = value
        if callable(self.value) and not dim:
            raise ValueError(
                "When using a flux callback the dimension has to be specified")
        elif not callable(self.value):
            dim = len(value)
        self.dim = dim
        self.data_type = data_type

    @property
    def value_is_callable(self):
        return callable(self.value)

    @property
    def additional_data(self):
        if self.value_is_callable:
            return [(f'flux_{i}', ps.typing.BasicType(
                self.data_type)) for i in range(self.dim)]
        else:
            return []

    @property
    def additional_data_init_callback(self):
        if self.value_is_callable:
            return self.value

    def __call__(self, field, direction_symbol, index_field, **kwargs):
        assert ps.FieldType.is_staggered(field)

        value = [index_field(f'flux_{i}') for i in range(
            self.dim)] if self.value_is_callable else self.value
        value = sp.Matrix(value)

        assert all([s == 0 for s in self.stencil[0]])
        accesses = [field.staggered_access(
            ps.stencil.offset_to_direction_string(d)) for d in self.stencil[1:]]

        conds = [
            sp.Equality(
                direction_symbol,
                ps.typing.CastFunc(
                    d + 1,
                    np.int32)) for d in range(
                len(accesses))]

        # use conditional
        conditional = None
        for access, condition, direction in zip(
                accesses, conds, self.stencil[1:]):
            d = sp.Matrix(direction)

            local_value = value

            # make sure the vector-access is non-negative
            if isinstance(access, sp.Mul):
                access *= -1
                local_value *= -1

            assignment = [
                ps.Assignment(
                    access,
                    local_value.dot(d) /
                    self.stencil.D ** 2)]

            # build stacked if-conditions for directions
            conditional = ps.astnodes.Conditional(
                condition, ps.astnodes.Block(assignment), conditional)

        return [conditional]

    def __hash__(self):
        return hash((Flux, self.stencil, self.value))

    def __eq__(self, other):
        return isinstance(
            other, Flux) and other.stencil == self.stencil and self.value == other.value


class DirichletAdditionalDataHandler(
        pystencils_walberla.additional_data_handler.AdditionalDataHandler):
    def __init__(self, stencil, boundary_object):
        assert isinstance(boundary_object, ps.boundaries.Dirichlet)
        self._boundary_object = boundary_object
        assert boundary_object.data_type in ("float32", "float64", "double")
        self.data_type = "float" if boundary_object.data_type == "float32" else "double"
        super().__init__(stencil=stencil)

    @property
    def constructor_arguments(self):
        return f", std::function<{self.data_type}(const Cell &, const shared_ptr<StructuredBlockForest>&, IBlock&)>& " \
               "dirichletCallback "

    @property
    def initialiser_list(self):
        return "elementInitaliser(dirichletCallback),"

    @property
    def additional_arguments_for_fill_function(self):
        return "blocks, "

    @property
    def additional_parameters_for_fill_function(self):
        return " const shared_ptr<StructuredBlockForest> &blocks, "

    def data_initialisation(self, _):
        init_list = [f"{self.data_type} InitialisatonAdditionalData = elementInitaliser(Cell(it.x(), it.y(), it.z()), "
                     "blocks, *block);", "element.value = InitialisatonAdditionalData;"]

        return "\n".join(init_list)

    @property
    def additional_member_variable(self):
        return f"std::function<{self.data_type}(const Cell &, const shared_ptr<StructuredBlockForest>&, IBlock&)> " \
               "elementInitaliser; "


class FluxAdditionalDataHandler(
        pystencils_walberla.additional_data_handler.AdditionalDataHandler):
    def __init__(self, stencil, boundary_object):
        self._boundary_object = boundary_object
        assert boundary_object.data_type in ("float32", "float64", "double")
        self.data_type = "float" if boundary_object.data_type == "float32" else "double"
        super().__init__(stencil=stencil)

    @property
    def constructor_arguments(self):
        return f", std::function<Vector3<{self.data_type}>(const Cell &, const shared_ptr<StructuredBlockForest>&, IBlock&)>& " \
               "fluxCallback "

    @property
    def initialiser_list(self):
        return "elementInitaliser(fluxCallback),"

    @property
    def additional_arguments_for_fill_function(self):
        return "blocks, "

    @property
    def additional_parameters_for_fill_function(self):
        return " const shared_ptr<StructuredBlockForest> &blocks, "

    def data_initialisation(self, direction):
        dirVec = self.stencil_info[direction][1]

        init_list = [
            f"Vector3<{self.data_type}> InitialisatonAdditionalData = elementInitaliser(Cell(it.x() + {dirVec[0]}, it.y() + {dirVec[1]}, it.z() + {dirVec[2]}), "
            "blocks, *block);", "element.flux_0 = InitialisatonAdditionalData[0];",
            "element.flux_1 = InitialisatonAdditionalData[1];"]
        if self._dim == 3:
            init_list.append(
                "element.flux_2 = InitialisatonAdditionalData[2];")

        return "\n".join(init_list)

    @property
    def additional_member_variable(self):
        return f"std::function<Vector3<{self.data_type}>(const Cell &, const shared_ptr<StructuredBlockForest>&, IBlock&)> " \
               "elementInitaliser; "


# this custom boundary generator is necessary because our boundary
# condition writes to several fields at once which is impossible with the
# shipped one
def generate_boundary(
        generation_context,
        stencil,
        class_name,
        dim: int,
        assignment,
        target=ps.enums.Target.CPU,
        data_type=None,
        cpu_openmp=None,
        namespace="pystencils",
        interface_mappings=(),
        generate_functor=True,
        **create_kernel_params,
):
    struct_name = "IndexInfo"

    config = pystencils_walberla.codegen.config_from_context(
        generation_context,
        target=target,
        data_type=data_type,
        cpu_openmp=cpu_openmp,
        **create_kernel_params,
    )
    create_kernel_params = config.__dict__
    del create_kernel_params["target"]
    del create_kernel_params["index_fields"]

    coordinate_names = ("x", "y", "z")[:dim]

    index_struct_dtype = np.dtype(
        [(name, np.int32) for name in coordinate_names], align=True
    )

    index_field = ps.Field(
        "indexVector",
        ps.FieldType.INDEXED,
        index_struct_dtype,
        layout=[0],
        shape=(
            ps.typing.TypedSymbol(
                "indexVectorSize", ps.typing.BasicType(np.int32)
            ),
            1,
        ),
        strides=(1, 1),
    )

    kernel_config = ps.CreateKernelConfig(
        index_fields=[index_field], target=target, **create_kernel_params
    )

    kernel = ps.kernelcreation.create_kernel(assignment, config=kernel_config)

    if isinstance(kernel, ps.astnodes.KernelFunction):
        kernel.function_name = f"boundary_{class_name}"
        selection_tree = pystencils_walberla.kernel_selection.KernelCallNode(
            kernel)
    elif isinstance(kernel, pystencils_walberla.kernel_selection.AbstractKernelSelectionNode):
        selection_tree = kernel
    else:
        raise ValueError(
            f"kernel_creation_function returned wrong type: {kernel.__class__}"
        )

    kernel_family = pystencils_walberla.kernel_selection.KernelFamily(
        selection_tree, class_name)
    interface_spec = pystencils_walberla.kernel_selection.HighLevelInterfaceSpec(
        kernel_family.kernel_selection_parameters, interface_mappings
    )

    additional_data_handler = pystencils_walberla.additional_data_handler.AdditionalDataHandler(
        stencil=stencil)

    context = {
        "kernel": kernel_family,
        "class_name": class_name,
        "interface_spec": interface_spec,
        "generate_functor": generate_functor,
        "StructName": struct_name,
        "StructDeclaration": pystencils_walberla.boundary.struct_from_numpy_dtype(struct_name, index_struct_dtype),
        "dim": dim,
        "target": target.name.lower(),
        "namespace": namespace,
        "inner_or_boundary": False,
        "single_link": False,
        "additional_data_handler": additional_data_handler,
    }

    env = jinja2.Environment(
        loader=jinja2.PackageLoader("pystencils_walberla"), undefined=jinja2.StrictUndefined
    )
    pystencils_walberla.jinja_filters.add_pystencils_filters_to_jinja_env(env)
    custom_env = jinja2.Environment(
        loader=jinja2.FileSystemLoader(pathlib.Path(__file__).parent), undefined=jinja2.StrictUndefined
    )
    pystencils_walberla.jinja_filters.add_pystencils_filters_to_jinja_env(
        custom_env)

    header = custom_env.get_template(
        "templates/Boundary.tmpl.h").render(**context)
    source = env.get_template("Boundary.tmpl.cpp").render(**context)

    source_extension = "cpp" if target == ps.enums.Target.CPU else "cu"
    generation_context.write_file(f"{class_name}.h", header)
    generation_context.write_file(f"{class_name}.{source_extension}", source)


def generate_kernel_selector(
        generation_context,
        class_name,
        namespace="pystencils",
        max_num_reactants=None,
        precision_suffix=None,
):
    """
    Generate helper functions to select a kernel with the appropriate
    floating-point precision and number of ek species for the currently
    active ek reaction and ek lattice.
    """

    context = {
        "namespace": namespace,
        "class_name": class_name,
        "precision_suffix": precision_suffix,
        "max_num_reactants": max_num_reactants,
    }

    custom_env = jinja2.Environment(
        loader=jinja2.FileSystemLoader(pathlib.Path(__file__).parent),
        undefined=jinja2.StrictUndefined
    )

    header = custom_env.get_template(
        "templates/ReactionKernelSelector.tmpl.h").render(**context)

    generation_context.write_file(f"{class_name}_all.h", header)
