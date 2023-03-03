#
# Copyright (C) 2020-2022 The ESPResSo project
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

import setuptools

import numpy as np
import sympy as sp
from sympy.tensor import IndexedBase

import pystencils as ps
import pystencils.astnodes

from pystencils.typing import BasicType, CastFunc, TypedSymbol
import pystencils_walberla
import pystencils_espresso

import lbmpy
import lbmpy.creationfunctions
import lbmpy.forcemodels
import lbmpy.stencils
import lbmpy.enums

import lbmpy_walberla
import lbmpy_espresso

import relaxation_rates
import code_generation_context
from pystencils.backends.cuda_backend import CudaBackend, CustomSympyPrinter


target = ps.Target.GPU

# Make sure we have the correct versions of the required dependencies
SpecifierSet = setuptools.version.pkg_resources.packaging.specifiers.SpecifierSet
for module, requirement in [(ps, '==1.1.1'), (lbmpy, '==1.1.1')]:
    assert SpecifierSet(requirement).contains(module.__version__), \
        f"{module.__name__} version {module.__version__} doesn't match requirement {requirement}"





# taken form walberla_lbm_generation.py with modificaitons
def field_and_symbol_substitute(expr, variable_prefix="lm.", variables_without_prefix=None):
    if variables_without_prefix is None:
        variables_without_prefix = []
    variables_without_prefix = [v.name if isinstance(v, sp.Symbol) else v for v in variables_without_prefix]
    substitutions = {}
    # check for member access
    if variable_prefix.endswith('.'):
        postfix = '_'
    else:
        postfix = ''
    for sym in expr.atoms(sp.Symbol):
        if isinstance(sym, ps.Field.Access):
            fa = sym
            prefix = "" if fa.field.name in variables_without_prefix else variable_prefix
            if prefix.endswith('.'):
                postfix2 = '_'
            else:
                postfix2 = ''
            if fa.field.index_dimensions == 0:
                substitutions[fa] = sp.Symbol(f"{prefix}{fa.field.name + postfix2}->get(x,y,z)")
            else:
                assert fa.field.index_dimensions == 1, "walberla supports only 0 or 1 index dimensions"
                substitutions[fa] = sp.Symbol(f"{prefix}{fa.field.name + postfix2}->get(x,y,z,{fa.index[0]})")
        else:
            if sym.name not in variables_without_prefix:
                substitutions[sym] = sp.Symbol(variable_prefix + sym.name + postfix)
    return expr.subs(substitutions)


def type_expr(eq, dtype):
    # manually cast 0.5 to dtype since this is somehow not done automatically
    eq = eq.subs(sp.Rational(1, 2), CastFunc(sp.Rational(1, 2), dtype))
    return eq.subs({s: TypedSymbol(s.name, dtype) for s in eq.atoms(sp.Symbol)})
def equations_to_code(equations, variable_prefix="", variables_without_prefix=None, dtype=None,
    backend = CudaBackend()):
    if dtype is None:
        dtype = BasicType("float64")

    if variables_without_prefix is None:
        variables_without_prefix = []
    if isinstance(equations, ps.AssignmentCollection):
        equations = equations.all_assignments

    variables_without_prefix = list(variables_without_prefix)
    result = []
    left_hand_side_names = [e.lhs.name for e in equations]
    for eq in equations:
        assignment = ps.astnodes.SympyAssignment(type_expr(eq.lhs, dtype=dtype),
                                     type_expr(field_and_symbol_substitute(eq.rhs, variable_prefix,
                                                                           variables_without_prefix
                                                                           + left_hand_side_names),
                                               dtype=dtype))
        result.append(backend(assignment))
    return "\n".join(result)

with code_generation_context.CodeGeneration() as ctx:
    ctx.double_accuracy = False
    ctx.cuda = True
    if target == ps.Target.GPU:
        target_prefix = "Cuda"
    else:
        target_prefix = ""

    params = {"target": target}

    # codegen configuration
    config = pystencils_espresso.generate_config(ctx, params)

    precision_prefix = pystencils_espresso.precision_prefix[ctx.double_accuracy]
    precision_suffix = pystencils_espresso.precision_suffix[ctx.double_accuracy]
    precision_rng = pystencils_espresso.precision_rng[ctx.double_accuracy]
    kT = sp.symbols('kT')
    stencil = lbmpy.stencils.LBStencil(lbmpy.enums.Stencil.D3Q19)
    fields = pystencils_espresso.generate_fields(config, stencil)

    force_field = fields['force']

    # LB Method definition
    method = lbmpy.creationfunctions.create_mrt_orthogonal(
        stencil=stencil,
        compressible=True,
        weighted=True,
        relaxation_rates=relaxation_rates.rr_getter,
        force_model=lbmpy.forcemodels.Schiller(force_field.center_vector)
    )
    cqc = method.conserved_quantity_computation
   
    default_dtype = np.float32 
    vel_symbols =cqc.velocity_symbols 
    vel_in_symbols = sp.symbols(f'vel_in_[:{method.stencil.D}]')
    rho_sym = sp.Symbol("rho")
    pdfs_sym = sp.symbols(f'f[:{method.stencil.Q}]')
    vel_arr_symbols = [IndexedBase(TypedSymbol('u', default_dtype), shape=(1,))[i] for i in range(len(vel_symbols))]
    momentum_density_symbols = sp.symbols(f'md_:{len(vel_symbols)}')
    
    eq_input_from_input_eqs = cqc.equilibrium_input_equations_from_init_values(sp.Symbol('rho_in'), vel_arr_symbols)

    constant_suffix = "f" 


    assignments = {

        'density_getters': equations_to_code(cqc.output_equations_from_pdfs(pdfs_sym, {"density": rho_sym}),
                                             variables_without_prefix=[e.name for e in pdfs_sym], dtype=default_dtype),
        'momentum_density_getter': equations_to_code( 
            cqc.output_equations_from_pdfs(pdfs_sym, {'density': rho_sym,
                                                                        'momentum_density': momentum_density_symbols}),
                                             variables_without_prefix=[e.name for e in pdfs_sym], dtype=default_dtype),
            "pdfs_setter": 
               equations_to_code(
                lbmpy.macroscopic_value_kernels.macroscopic_values_setter(
        method,
        sp.Symbol("rho_in"),
        vel_in_symbols,
        pdfs_sym),
                                             variables_without_prefix=[e.name for e in pdfs_sym], dtype=default_dtype)
       }
        
 
    for k,v in assignments.items():
        print(f"{k}: \n{v}\n")

