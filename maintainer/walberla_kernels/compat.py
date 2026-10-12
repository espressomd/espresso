#
# Copyright (C) 2026 The ESPResSo project
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

"""
Generate sweepgen sweeps behind the C++ interface of the classes that
``pystencils_walberla.generate_sweep()`` used to produce, so that the
waLBerla bridge can use them unchanged.
"""

import re

import sweepgen
from pystencils.types import PsCustomType, PsPointerType
from pystencilssfg.lang import Ref

LEGACY_NAMESPACE = "walberla::pystencils"
INTERNAL_NAMESPACE = "internal"
FIELD_ID_TYPE = "walberla::BlockDataID"
BLOCK_OFFSET_PROPERTY = "block_offset"


def legacy_accessor_suffix(name):
    """Suffix of the legacy getter/setter of a parameter, e.g. ``kT`` -> ``Kt``."""
    return name.capitalize()


def _var(sfg, name, cpp_type):
    return sfg.var(name, PsCustomType(cpp_type))


def _cref(sfg, name, cpp_type):
    return sfg.var(name, Ref(PsCustomType(cpp_type, const=True)))


def legacy_sweep(sfg, ctx, name, assignments, *, ctor_order=None,
                 gen_config=None, swap_fields=(), legacy_types=None,
                 athermal_assignments=None, athermal_param="kT"):
    """
    Generate the sweep ``name`` as a sweepgen class in namespace
    ``walberla::pystencils::internal`` and a wrapper class
    ``walberla::pystencils::<name>`` with the legacy interface.

    The wrapper constructor takes the field IDs in alphabetical order,
    then the scalar parameters in alphabetical order, unless
    ``ctor_order`` lists the parameter names explicitly.

    ``swap_fields`` is a list of ``(field, shadow_field)`` pairs that the
    sweep swaps after each run (``Sweep.swap_fields()``).

    A sweep property ``block_offset`` (``Vector3< int64_t >``, e.g. the
    offset of a counter-based RNG) isn't a constructor argument of the
    wrapper: it starts at zero and ``configure()`` sets it to the lower
    corner of the block cell bounding box. It is exposed as
    ``get/setBlock_offset_<i>``.

    ``legacy_types`` maps scalar names to the C++ integer type of the
    wrapper API when it differs from the type in the sweep (e.g.
    ``{"time_step": "uint32_t"}``); the wrapper converts on access.

    If ``athermal_assignments`` is given, a second sweep
    ``<name>Athermal`` is generated from it (same fields and swaps), and
    ``run()`` dispatches to it when the parameter ``athermal_param`` of the
    first sweep is not positive. This avoids evaluating an RNG in
    deterministic runs. Scalars of both sweeps are set together.
    """
    sweep = sweepgen.Sweep(name, assignments, gen_config)
    if sweep.sparse:
        raise NotImplementedError(
            f"{name}: sparse sweeps are not supported by the legacy wrapper "
            f"yet")
    if athermal_assignments is not None and not swap_fields:
        raise NotImplementedError(
            f"{name}: athermal dispatch needs a swapping sweep")
    for field, shadow_field in swap_fields:
        sweep.swap_fields(field, shadow_field)
    athermal = None
    if athermal_assignments is not None:
        athermal = sweepgen.Sweep(f"{name}Athermal", athermal_assignments,
                                  gen_config)
        for field, shadow_field in swap_fields:
            athermal.swap_fields(field, shadow_field)
    with sfg.namespace(LEGACY_NAMESPACE):
        with sfg.namespace(INTERNAL_NAMESPACE):
            sfg.generate(sweep)
            if athermal is not None:
                sfg.generate(athermal)
        _legacy_wrapper(sfg, ctx, name, sweep, ctor_order,
                        legacy_types or {}, bool(swap_fields), athermal,
                        athermal_param)
    return sweep


def _legacy_wrapper(sfg, ctx, name, sweep, ctor_order, legacy_types,
                    has_swaps, athermal=None, athermal_param="kT"):
    container = sweep.properties_container
    inner_params = list(container.ctor_params)
    properties = container.properties
    offset_prop = properties.get(BLOCK_OFFSET_PROPERTY)
    offset_type = None
    if offset_prop is not None:
        offset_type = offset_prop.dtype.c_string()
        match = re.fullmatch(r"walberla::Vector3< (u?int(?:32|64)_t) >",
                             offset_type)
        if match is None:
            raise NotImplementedError(
                f"{name}: property {BLOCK_OFFSET_PROPERTY} of type "
                f"{offset_type} is not supported by the legacy wrapper yet")
        offset_elem = match.group(1)
    field_ids = sorted(
        p.name for p in properties.values()
        if not p.setter and p.dtype.c_string() == FIELD_ID_TYPE)
    scalars = sorted(p.name for p in properties.values()
                     if p.setter and p.name != BLOCK_OFFSET_PROPERTY)
    supported = set(field_ids) | set(scalars)
    if offset_prop is not None:
        supported.add(BLOCK_OFFSET_PROPERTY)
    unsupported = ({p.name for p in properties.values()} |
                   {p.name for p in inner_params}) - supported
    if unsupported:
        raise NotImplementedError(
            f"{name}: sweep properties {sorted(unsupported)} are not "
            f"supported by the legacy wrapper yet")
    for key, value in legacy_types.items():
        if key not in scalars and not (
                key == BLOCK_OFFSET_PROPERTY and offset_prop is not None):
            raise ValueError(f"{name}: legacy type for unknown scalar {key}")
        if re.fullmatch(r"u?int(8|16|32|64)_t", value) is None:
            raise ValueError(
                f"{name}: unsupported legacy type {value} for {key}")
    order = list(ctor_order) if ctor_order is not None else field_ids + scalars
    ctor_names = [p.name for p in inner_params
                  if p.name != BLOCK_OFFSET_PROPERTY]
    if sorted(order) != sorted(ctor_names):
        raise ValueError(
            f"constructor order {order} of {name} doesn't match the sweep "
            f"parameters {ctor_names}")
    params_by_name = {p.name: p for p in inner_params}

    def api_type(scalar):
        return legacy_types.get(scalar, properties[scalar].dtype.c_string())

    # wrapper constructor parameters (with legacy types) and the arguments
    # they are forwarded as to the sweepgen sweep
    wrapper_params = {
        n: (_var(sfg, n, legacy_types[n]) if n in legacy_types
            else params_by_name[n]) for n in order}
    forwarded = [
        f"{offset_type}(0)" if p.name == BLOCK_OFFSET_PROPERTY
        else wrapper_params[p.name] for p in inner_params]

    impl = _var(sfg, "impl_", f"{INTERNAL_NAMESPACE}::{name}")
    impl_athermal = None
    athermal_scalars = set()
    if athermal is not None:
        if athermal_param not in scalars:
            raise ValueError(
                f"{name}: athermal parameter {athermal_param} is not a "
                f"scalar of the sweep")
        impl_athermal = _var(sfg, "impl_athermal_",
                             f"{INTERNAL_NAMESPACE}::{name}Athermal")
        athermal_container = athermal.properties_container
        athermal_params = []
        for p in athermal_container.ctor_params:
            if p.name not in wrapper_params:
                raise ValueError(
                    f"{name}: athermal sweep parameter {p.name} has no "
                    f"counterpart in the sweep")
            athermal_params.append(wrapper_params[p.name])
        athermal_scalars = {p.name for p in athermal_container.properties.values()
                            if p.setter}
        if not athermal_scalars <= set(scalars):
            raise ValueError(f"{name}: athermal sweep has scalars "
                             f"{sorted(athermal_scalars - set(scalars))} "
                             f"that the sweep doesn't have")
    block = sfg.var("block", PsPointerType(PsCustomType("walberla::IBlock")))
    blocks = _cref(sfg, "blocks",
                   "walberla::shared_ptr< walberla::StructuredBlockStorage >")
    ci = _cref(sfg, "globalCellInterval", "walberla::CellInterval")
    ghost_layers = _var(sfg, "ghostLayers", "walberla::cell_idx_t")
    kernel = _cref(sfg, "kernel", f"walberla::shared_ptr< {name} >")

    sfg.include("<functional>")
    sfg.include("core/DataTypes.h")
    sfg.include("core/cell/CellInterval.h")
    sfg.include("domain_decomposition/IBlock.h")
    sfg.include("domain_decomposition/StructuredBlockStorage.h")

    if ctx.is_gpu:
        sfg.include("gpu/GPUWrapper.h")
        stream = _var(sfg, "stream", "gpuStream_t")
        stream_param = [(stream, "nullptr")]
        exec_tag = ", walberla::v8::sweep::exectag::GPU{stream}"
        ci_sweep_fn = "std::function< void(walberla::IBlock *, gpuStream_t) >"
        run_args = "block, stream"
        lambda_stream = ", gpuStream_t stream = nullptr"
        capture_stream = ", stream"
    else:
        stream_param = []
        exec_tag = ""
        ci_sweep_fn = "std::function< void(walberla::IBlock *) >"
        run_args = "block"
        lambda_stream = ""
        capture_stream = ""
    sweep_fn = "std::function< void(walberla::IBlock *) >"

    if athermal is None:
        run_body = (f"impl_(block{exec_tag});",)
    else:
        run_body = (f"if (impl_.{athermal_param}() > 0) {{",
                    f"  impl_(block{exec_tag});",
                    "} else {",
                    f"  impl_athermal_(block{exec_tag});",
                    "}")
    methods = [
        sfg.method("run").params(block, *stream_param).inline()(*run_body),
        sfg.method("operator()").params(block, *stream_param).inline()(
            f"run({run_args});"),
        sfg.method("getSweep").params(kernel).returns(
            PsCustomType(sweep_fn)).static().inline()(
            "return [kernel](walberla::IBlock *b) { kernel->run(b); };"),
        sfg.method("getSweep").params(*stream_param).returns(
            PsCustomType(sweep_fn)).inline()(
            f"return [this{capture_stream}](walberla::IBlock *b) "
            f"{{ this->run(b{capture_stream}); }};"),
    ]
    if not has_swaps:
        # runOnCellInterval() doesn't swap the fields of a swapping sweep
        # (the shadow field is private to the sweep), so it is only
        # provided for sweeps without swaps
        methods += [
            sfg.method("runOnCellInterval").params(
                blocks, ci, ghost_layers, block, *stream_param).inline()(
                "walberla::CellInterval ci = globalCellInterval;",
                "walberla::CellInterval blockBB = "
                "blocks->getBlockCellBB(*block);",
                "blockBB.expand(ghostLayers);",
                "ci.intersect(blockBB);",
                "blocks->transformGlobalToBlockLocalCellInterval(ci, *block);",
                "if (ci.empty())",
                "  return;",
                f"impl_.runOnCellInterval(block, ci{exec_tag});"),
            sfg.method("getSweepOnCellInterval").params(
                kernel, blocks, ci, (ghost_layers, "1")).returns(
                PsCustomType(ci_sweep_fn)).static().inline()(
                f"return [kernel, blocks, globalCellInterval, ghostLayers]"
                f"(walberla::IBlock *b{lambda_stream}) {{",
                f"  kernel->runOnCellInterval(blocks, globalCellInterval, "
                f"ghostLayers, b{capture_stream});",
                "};"),
            sfg.method("getSweepOnCellInterval").params(
                blocks, ci, (ghost_layers, "1"), *stream_param).returns(
                PsCustomType(sweep_fn)).inline()(
                f"return [this, blocks, globalCellInterval, ghostLayers"
                f"{capture_stream}](walberla::IBlock *b) {{",
                f"  this->runOnCellInterval(blocks, globalCellInterval, "
                f"ghostLayers, b{capture_stream});",
                "};"),
        ]
    if offset_prop is None:
        methods.append(sfg.method("configure").params(
            _cref(sfg, "", "walberla::shared_ptr< walberla::StructuredBlockStorage >"),
            sfg.var("", PsPointerType(PsCustomType("walberla::IBlock")))
        ).inline()())
    else:
        methods.append(sfg.method("configure").params(blocks, block).inline()(
            "const walberla::Cell blockMin = "
            "blocks->getBlockCellBB(*block).min();",
            f"impl_.{BLOCK_OFFSET_PROPERTY}() = "
            f"{offset_type}(blockMin[0], blockMin[1], blockMin[2]);"))
        offset_api = legacy_types.get(BLOCK_OFFSET_PROPERTY, offset_elem)
        for i in range(3):
            suffix = legacy_accessor_suffix(f"{BLOCK_OFFSET_PROPERTY}_{i}")
            value = sfg.var("value", PsCustomType(offset_api))
            methods.append(
                sfg.method(f"get{suffix}").returns(
                    PsCustomType(offset_api)).const().inline()(
                    f"return static_cast< {offset_api} >"
                    f"(impl_.{BLOCK_OFFSET_PROPERTY}()[{i}]);"))
            methods.append(
                sfg.method(f"set{suffix}").params(value).inline()(
                    f"impl_.{BLOCK_OFFSET_PROPERTY}()[{i}] = value;"))
    for scalar in scalars:
        dtype = properties[scalar].dtype
        cpp_type = api_type(scalar)
        api_dtype = PsCustomType(cpp_type)
        suffix = legacy_accessor_suffix(scalar)
        value = sfg.var("value", api_dtype)
        getter = (f"return static_cast< {cpp_type} >(impl_.{scalar}());"
                  if scalar in legacy_types else f"return impl_.{scalar}();")
        methods.append(
            sfg.method(f"get{suffix}").returns(api_dtype).const().inline()(
                getter))
        setter_body = [f"impl_.{scalar}() = value;"]
        if scalar in athermal_scalars:
            setter_body.append(f"impl_athermal_.{scalar}() = value;")
        methods.append(
            sfg.method(f"set{suffix}").params(value).inline()(*setter_body))

    constructor = sfg.constructor(*(wrapper_params[n] for n in order))
    constructor = constructor.init(impl)(*forwarded)
    members = [impl]
    if athermal is not None:
        constructor = constructor.init(impl_athermal)(*athermal_params)
        members.append(impl_athermal)
    sfg.klass(name)(
        sfg.public(constructor, *methods),
        sfg.private(*members),
    )
