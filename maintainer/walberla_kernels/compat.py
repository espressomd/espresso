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

import sweepgen
from pystencils.types import PsCustomType, PsPointerType
from pystencilssfg.lang import Ref

LEGACY_NAMESPACE = "walberla::pystencils"
INTERNAL_NAMESPACE = "internal"
FIELD_ID_TYPE = "walberla::BlockDataID"


def legacy_accessor_suffix(name):
    """Suffix of the legacy getter/setter of a parameter, e.g. ``kT`` -> ``Kt``."""
    return name.capitalize()


def _var(sfg, name, cpp_type):
    return sfg.var(name, PsCustomType(cpp_type))


def _cref(sfg, name, cpp_type):
    return sfg.var(name, Ref(PsCustomType(cpp_type, const=True)))


def legacy_sweep(sfg, ctx, name, assignments, *, ctor_order=None,
                 gen_config=None):
    """
    Generate the sweep ``name`` as a sweepgen class in namespace
    ``walberla::pystencils::internal`` and a wrapper class
    ``walberla::pystencils::<name>`` with the legacy interface.

    The wrapper constructor takes the field IDs in alphabetical order,
    then the scalar parameters in alphabetical order, unless
    ``ctor_order`` lists the parameter names explicitly.
    """
    sweep = sweepgen.Sweep(name, assignments, gen_config)
    if sweep.sparse:
        raise NotImplementedError(
            f"{name}: sparse sweeps are not supported by the legacy wrapper "
            f"yet")
    with sfg.namespace(LEGACY_NAMESPACE):
        with sfg.namespace(INTERNAL_NAMESPACE):
            sfg.generate(sweep)
        _legacy_wrapper(sfg, ctx, name, sweep, ctor_order)
    return sweep


def _legacy_wrapper(sfg, ctx, name, sweep, ctor_order):
    container = sweep.properties_container
    inner_params = list(container.ctor_params)
    field_ids = sorted(
        p.name for p in container.properties.values()
        if not p.setter and p.dtype.c_string() == FIELD_ID_TYPE)
    scalars = sorted(p.name for p in container.properties.values()
                     if p.setter)
    supported = set(field_ids) | set(scalars)
    unsupported = ({p.name for p in container.properties.values()} |
                   {p.name for p in inner_params}) - supported
    if unsupported:
        raise NotImplementedError(
            f"{name}: sweep properties {sorted(unsupported)} are not "
            f"supported by the legacy wrapper yet")
    order = list(ctor_order) if ctor_order is not None else field_ids + scalars
    if sorted(order) != sorted(p.name for p in inner_params):
        raise ValueError(
            f"constructor order {order} of {name} doesn't match the sweep "
            f"parameters {[p.name for p in inner_params]}")
    params_by_name = {p.name: p for p in inner_params}

    impl = _var(sfg, "impl_", f"{INTERNAL_NAMESPACE}::{name}")
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

    methods = [
        sfg.method("run").params(block, *stream_param).inline()(
            f"impl_(block{exec_tag});"),
        sfg.method("runOnCellInterval").params(
            blocks, ci, ghost_layers, block, *stream_param).inline()(
            "walberla::CellInterval ci = globalCellInterval;",
            "walberla::CellInterval blockBB = blocks->getBlockCellBB(*block);",
            "blockBB.expand(ghostLayers);",
            "ci.intersect(blockBB);",
            "blocks->transformGlobalToBlockLocalCellInterval(ci, *block);",
            "if (ci.empty())",
            "  return;",
            f"impl_.runOnCellInterval(block, ci{exec_tag});"),
        sfg.method("operator()").params(block, *stream_param).inline()(
            f"run({run_args});"),
        sfg.method("getSweep").params(kernel).returns(
            PsCustomType(sweep_fn)).static().inline()(
            "return [kernel](walberla::IBlock *b) { kernel->run(b); };"),
        sfg.method("getSweepOnCellInterval").params(
            kernel, blocks, ci, (ghost_layers, "1")).returns(
            PsCustomType(ci_sweep_fn)).static().inline()(
            f"return [kernel, blocks, globalCellInterval, ghostLayers]"
            f"(walberla::IBlock *b{lambda_stream}) {{",
            f"  kernel->runOnCellInterval(blocks, globalCellInterval, "
            f"ghostLayers, b{capture_stream});",
            "};"),
        sfg.method("getSweep").params(*stream_param).returns(
            PsCustomType(sweep_fn)).inline()(
            f"return [this{capture_stream}](walberla::IBlock *b) "
            f"{{ this->run(b{capture_stream}); }};"),
        sfg.method("getSweepOnCellInterval").params(
            blocks, ci, (ghost_layers, "1"), *stream_param).returns(
            PsCustomType(sweep_fn)).inline()(
            f"return [this, blocks, globalCellInterval, ghostLayers"
            f"{capture_stream}](walberla::IBlock *b) {{",
            f"  this->runOnCellInterval(blocks, globalCellInterval, "
            f"ghostLayers, b{capture_stream});",
            "};"),
        sfg.method("configure").params(
            _cref(sfg, "", "walberla::shared_ptr< walberla::StructuredBlockStorage >"),
            sfg.var("", PsPointerType(PsCustomType("walberla::IBlock")))
        ).inline()(),
    ]
    for scalar in scalars:
        dtype = container.properties[scalar].dtype
        suffix = legacy_accessor_suffix(scalar)
        value = sfg.var("value", dtype)
        methods.append(
            sfg.method(f"get{suffix}").returns(dtype).const().inline()(
                f"return impl_.{scalar}();"))
        methods.append(
            sfg.method(f"set{suffix}").params(value).inline()(
                f"impl_.{scalar}() = value;"))

    sfg.klass(name)(
        sfg.public(
            sfg.constructor(*(params_by_name[n] for n in order))
            .init(impl)(*inner_params),
            *methods,
        ),
        sfg.private(impl),
    )
