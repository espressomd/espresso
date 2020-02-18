# Copyright (C) 2019 The ESPResSo project
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

import re
import os
import io
import sys
import ast
import tokenize
import unittest
import importlib
import espressomd
from unittest.mock import MagicMock


def _id(x):
    return x


# global variable: if one import failed, all subsequent imports will be skipped,
# see skip_future_imports_dependency()
skip_future_imports = False


def configure_and_import(filepath,
                         gpu=False,
                         substitutions=lambda x: x,
                         cmd_arguments=None,
                         script_suffix=None,
                         move_to_script_dir=True,
                         random_seeds=True,
                         mock_visualizers=True,
                         **parameters):
    """
    Copy a Python script to a new location and alter some lines of code:

    - change global variables and local variables (up to 1 indentation level)
    - pass command line arguments during import to emulate shell execution
    - disable the OpenGL/Mayavi modules if they are not compiled
    - disable the matplotlib GUI using a text-based backend
    - use random seeds for the RNG in NumPy
    - temporarily move to the directory where the script is located

    Parameters
    ----------
    filepath : str
        python script to import
    gpu : bool
        whether GPU is necessary or not
    substitutions : function
        custom text replacement operation (useful to edit out calls to the
        OpenGL or Mayavi visualizers' ``run()`` method)
    cmd_arguments : list
        command line arguments, i.e. sys.argv without the script path
    script_suffix : str
        suffix to append to the configured script (useful when a single
        module is being tested by multiple tests in parallel)
    random_seeds : bool
        if ``True``, use random seeds in RNGs
    mock_visualizers : bool
        if ``True``, substitute ES visualizers with `Mock()` classes in case
        of `ImportError()` (use ``False`` if an `ImportError()` is relevant
        to your test)
    move_to_script_dir : bool
        if ``True``, move to the script's directory (useful when the script
        needs to load files hardcoded as relative paths, or when files are
        generated and need cleanup); this is enabled by default
    \*\*parameters :
        global variables to replace

    """
    if skip_future_imports:
        module = MagicMock()
        skipIfMissingImport = skip_future_imports_dependency(filepath)
        return module, skipIfMissingImport
    if gpu and not espressomd.gpu_available():
        skip_future_imports_dependency(filepath)
        skipIfMissingGPU = unittest.skip("gpu not available, skipping test!")
        module = MagicMock()
        return module, skipIfMissingGPU
    filepath = os.path.abspath(filepath)
    # load original script
    # read in binary mode, then decode as UTF-8 to avoid this python3.5 error:
    # UnicodeDecodeError: 'ascii' codec can't decode byte 0xc3 in position 915:
    #                      ordinal not in range(128)
    with open(filepath, "rb") as f:
        code = f.read().decode(encoding="utf-8")
    # custom substitutions
    code = substitutions(code)
    assert code.strip()
    # substitute global variables
    code = substitute_variable_values(code, **parameters)
    # substitute command line arguments
    if cmd_arguments is not None:
        code, old_sys_argv = set_cmd(code, filepath, cmd_arguments)
    # disable matplotlib GUI using the Agg backend
    code = disable_matplotlib_gui(code)
    # disable OpenGL/Mayavi GUI using MagicMock()
    if mock_visualizers:
        code = mock_es_visualization(code)
    # use random seeds for NumPy RNG
    if random_seeds:
        code = set_random_seeds(code)
    # save changes to a new file
    if script_suffix:
        if script_suffix[0] != "_":
            script_suffix = "_" + script_suffix
    else:
        script_suffix = ""
    script_suffix += "_processed.py"
    output_filepath = os.path.splitext(filepath)[0] + script_suffix
    assert os.path.isfile(output_filepath) is False, \
        "File {} already processed, cannot overwrite".format(output_filepath)
    with open(output_filepath, "wb") as f:
        f.write(code.encode(encoding="utf-8"))
    # import
    dirname, basename = os.path.split(output_filepath)
    if move_to_script_dir:
        os.chdir(dirname)
    sys.path.insert(0, dirname)
    module_name = os.path.splitext(basename)[0]
    try:
        module = importlib.import_module(module_name)
    except espressomd.FeaturesError as err:
        skip_future_imports_dependency(filepath)
        skipIfMissingFeatures = unittest.skip(str(err) + ", skipping test!")
        module = MagicMock()
    else:
        skipIfMissingFeatures = _id
    if cmd_arguments is not None:
        # restore original command line arguments
        sys.argv = old_sys_argv
    return module, skipIfMissingFeatures


class GetSysArgparseImports(ast.NodeVisitor):
    """
    Find all line numbers where ``sys`` or ``argparse`` is imported.
    """

    def __init__(self):
        self.linenos = []

    def visit_Import(self, node):
        for child in node.names:
            if child.name.split(".")[0] in ("sys", "argparse"):
                self.linenos.append(node.lineno)

    def visit_ImportFrom(self, node):
        if node.module.split(".")[0] in ("sys", "argparse"):
            self.linenos.append(node.lineno)


def set_cmd(code, filepath, cmd_arguments):
    """
    Set the ``sys.argv`` in the imported module while checkpointing
    the current ``sys.argv`` value.
    """
    assert isinstance(cmd_arguments, (list, tuple))
    sys_argv = list(map(str, cmd_arguments))
    sys_argv.insert(0, os.path.basename(filepath))
    visitor = GetSysArgparseImports()
    visitor.visit(ast.parse(protect_ipython_magics(code)))
    assert visitor.linenos, "module sys (or argparse) is not imported"
    lines = code.split("\n")
    mapping = delimit_statements(code)
    lineno = mapping[min(visitor.linenos)]
    line = lines[lineno - 1]
    indentation = line[:len(line) - len(line.lstrip())]
    lines[lineno - 1] = indentation + "import sys;sys.argv = " + \
        str(sys_argv) + ";" + line.lstrip()
    code = "\n".join(lines)
    old_sys_argv = list(sys.argv)
    return code, old_sys_argv


class GetVariableAssignments(ast.NodeVisitor):
    """
    Find all assignments in the global namespace.
    """

    def __init__(self, variable_names):
        self.variables = {x: [] for x in variable_names}

    def visit_Assign(self, node):
        for target in node.targets:
            if hasattr(target, "id"):
                varname = target.id
                if varname in self.variables:
                    assert len(node.targets) == 1, \
                        "Cannot substitute multiple assignments (variable " \
                        "'{}' at line {})".format(varname, node.lineno)
                    self.variables[varname].append(node.lineno)

    def visit_ClassDef(self, node):
        pass  # skip class scopes

    def visit_FunctionDef(self, node):
        pass  # skip function scopes

    def visit_AsyncFunctionDef(self, node):
        pass  # skip function scopes


def substitute_variable_values(code, strings_as_is=False, keep_original=True,
                               **parameters):
    """
    Substitute values of global variables.

    Parameters
    ----------
    code : str
        Source code to edit.
    strings_as_is : bool
        If ``True``, consider all values in \*\*parameters are strings and
        substitute them in-place without formatting by ``repr()``.
    keep_original : bool
        Keep the original value (e.g. ``N = 10; _N__original = 1000``), helps
        with debugging.
    \*\*parameters :
        Variable names and their new value.

    """
    tree = ast.parse(protect_ipython_magics(code))
    visitor = GetVariableAssignments(parameters.keys())
    visitor.visit(tree)
    # split lines
    lines = code.split("\n")
    mapping = delimit_statements(code)
    # substitute values
    for varname, new_value in parameters.items():
        linenos = visitor.variables[varname]
        assert linenos, "variable {} has no assignment".format(varname)
        new_value = strings_as_is and new_value or repr(new_value)
        for lineno in linenos:
            identation, old_value = lines[lineno - 1].split(varname, 1)
            lines[lineno - 1] = identation + varname + " = " + new_value
            if keep_original:
                lines[lineno - 1] += "; _" + varname + "__original" + old_value
            else:
                for lineno in range(lineno + 1, mapping[lineno]):
                    lines[lineno - 1] = ""
    return "\n".join(lines)


class GetPrngSeedEspressomdSystem(ast.NodeVisitor):
    """
    Find all assignments of :class:`espressomd.system.System` in the global
    namespace. Assignments made in classes or function raise an error. Detect
    random seed setup of the numpy PRNG.
    """

    def __init__(self):
        self.numpy_random_aliases = set()
        self.es_system_aliases = set()
        self.variable_system_aliases = set()
        self.numpy_seeds = []
        self.abort_message = None
        self.error_msg_multi_assign = "Cannot parse {} in a multiple assignment (line {})"

    def visit_Import(self, node):
        # get system aliases
        for child in node.names:
            if child.name == "espressomd.system.System":
                name = (child.asname or child.name)
                self.es_system_aliases.add(name)
            elif child.name == "espressomd.system":
                name = (child.asname or child.name)
                self.es_system_aliases.add(name + ".System")
            elif child.name == "espressomd.System":
                name = (child.asname or child.name)
                self.es_system_aliases.add(name)
            elif child.name == "espressomd":
                name = (child.asname or "espressomd")
                self.es_system_aliases.add(name + ".system.System")
                self.es_system_aliases.add(name + ".System")
            elif child.name == "numpy.random":
                name = (child.asname or child.name)
                self.numpy_random_aliases.add(name)
            elif child.name == "numpy":
                name = (child.asname or "numpy")
                self.numpy_random_aliases.add(name + ".random")

    def visit_ImportFrom(self, node):
        # get system aliases
        for child in node.names:
            if node.module == "espressomd" and child.name == "system":
                name = (child.asname or child.name)
                self.es_system_aliases.add(name + ".System")
            elif node.module == "espressomd" and child.name == "System":
                self.es_system_aliases.add(child.asname or child.name)
            elif node.module == "espressomd.system" and child.name == "System":
                self.es_system_aliases.add(child.asname or child.name)
            elif node.module == "numpy" and child.name == "random":
                self.numpy_random_aliases.add(child.asname or child.name)
            elif node.module == "numpy.random":
                self.numpy_random_aliases.add(child.asname or child.name)

    def is_es_system(self, child):
        if hasattr(child, "value"):
            if hasattr(child.value, "value") and hasattr(
                    child.value.value, "id"):
                if (child.value.value.id + "." + child.value.attr +
                        "." + child.attr) in self.es_system_aliases:
                    return True
            else:
                if hasattr(child.value, "id") and (child.value.id + "." +
                                                   child.attr) in self.es_system_aliases:
                    return True
        elif isinstance(child, ast.Call):
            if hasattr(child, "id") and child.func.id in self.es_system_aliases:
                return True
            elif hasattr(child, "func") and hasattr(child.func, "value") and (
                    hasattr(child.func.value, "value") and
                    (child.func.value.value.id + "." +
                     child.func.value.attr + "." + child.func.attr)
                    or (child.func.value.id + "." + child.func.attr)) in self.es_system_aliases:
                return True
            elif hasattr(child, "func") and hasattr(child.func, "id") and child.func.id in self.es_system_aliases:
                return True
        elif hasattr(child, "id") and child.id in self.es_system_aliases:
            return True
        return False

    def detect_es_system_instances(self, node):
        varname = None
        for target in node.targets:
            if isinstance(target, ast.Name):
                if hasattr(target, "id") and hasattr(node.value, "func") and \
                        self.is_es_system(node.value.func):
                    varname = target.id
            elif isinstance(target, ast.Tuple):
                value = node.value
                if (isinstance(value, ast.Tuple) or isinstance(value, ast.List)) \
                        and any(map(self.is_es_system, node.value.elts)):
                    raise AssertionError(self.error_msg_multi_assign.format(
                        "espressomd.System", node.lineno))
        if varname is not None:
            assert len(node.targets) == 1, self.error_msg_multi_assign.format(
                "espressomd.System", node.lineno)
            assert self.abort_message is None, \
                "Cannot process espressomd.System assignments in " + self.abort_message
            self.variable_system_aliases.add(varname)

    def detect_np_random_expr_seed(self, node):
        if hasattr(node.value, "func") and hasattr(node.value.func, "value") \
                and (hasattr(node.value.func.value, "id") and node.value.func.value.id in self.numpy_random_aliases
                     or hasattr(node.value.func.value, "value")
                     and hasattr(node.value.func.value.value, "id")
                     and hasattr(node.value.func.value, "attr")
                     and (node.value.func.value.value.id + "." +
                          node.value.func.value.attr) in self.numpy_random_aliases
                     ) and node.value.func.attr == "seed":
            self.numpy_seeds.append(node.lineno)

    def visit_Assign(self, node):
        self.detect_es_system_instances(node)

    def visit_Expr(self, node):
        self.detect_np_random_expr_seed(node)

    def visit_ClassDef(self, node):
        self.abort_message = "class definitions"
        ast.NodeVisitor.generic_visit(self, node)
        self.abort_message = None

    def visit_FunctionDef(self, node):
        self.abort_message = "function definitions"
        ast.NodeVisitor.generic_visit(self, node)
        self.abort_message = None

    def visit_AsyncFunctionDef(self, node):
        self.abort_message = "function definitions"
        ast.NodeVisitor.generic_visit(self, node)
        self.abort_message = None


def set_random_seeds(code):
    """
    Remove numpy random number generator seeds, such that
    the sample/tutorial always starts with different random seeds.
    """
    code = protect_ipython_magics(code)
    tree = ast.parse(code)
    visitor = GetPrngSeedEspressomdSystem()
    visitor.visit(tree)
    lines = code.split("\n")
    # delete explicit NumPy seed
    for lineno in visitor.numpy_seeds:
        old_stmt = ".seed"
        new_stmt = ".seed;_random_seed_np = (lambda *args, **kwargs: None)"
        lines[lineno - 1] = lines[lineno - 1].replace(old_stmt, new_stmt, 1)
    code = "\n".join(lines)
    code = deprotect_ipython_magics(code)
    return code


def delimit_statements(code):
    """
    For every Python statement, map the line number where it starts to the
    line number where it ends.
    """
    statements = []
    statement_start = None
    for tok in tokenize.tokenize(io.BytesIO(code.encode("utf-8")).readline):
        if tok.exact_type == tokenize.ENDMARKER:
            break
        elif tok.exact_type == tokenize.ENCODING:
            pass
        elif tok.start == tok.end:
            pass
        elif tok.exact_type == tokenize.NEWLINE or tok.exact_type == tokenize.NL and prev_tok.exact_type == tokenize.COMMENT:
            statements.append((statement_start, tok.start[0]))
            statement_start = None
        elif tok.exact_type == tokenize.NL:
            pass
        elif statement_start is None:
            statement_start = tok.start[0]
        prev_tok = tok
    return dict(statements)


def protect_ipython_magics(code):
    """
    Replace all IPython magic commands (e.g. ``%matplotlib notebook``) by a
    formatted comment. This is necessary whenever the code must be parsed
    through ``ast``, because magic commands are not valid Python statements.
    """
    return re.sub("^(%+)(?=[a-z])", "#_IPYTHON_MAGIC_\g<1>", code, flags=re.M)


def deprotect_ipython_magics(code):
    """
    Reverse the action of :func:`protect_ipython_magics`.
    """
    return re.sub("^#_IPYTHON_MAGIC_(%+)(?=[a-z])", "\g<1>", code, flags=re.M)


class GetMatplotlibPyplot(ast.NodeVisitor):
    """
    Find all line numbers where ``matplotlib`` and ``pyplot`` are imported
    and store their alias. Find line numbers where the matplotlib backend
    is set, where the pyplot interactive mode is activated and where IPython
    magic commands (in the processed form ``get_ipython().func(args)``) are
    set.
    """

    def __init__(self):
        self.matplotlib_first = None
        self.matplotlib_aliases = []
        self.pyplot_aliases = []
        self.pyplot_paths = set()
        self.pyplot_interactive_linenos = []
        self.matplotlib_backend_linenos = []
        self.ipython_magic_linenos = []

    def make_pyplot_paths(self):
        self.pyplot_paths = set(tuple(x.split("."))
                                for x in self.pyplot_aliases)

    def visit_Import(self, node):
        # get line number of the first matplotlib import
        for child in node.names:
            if child.name.split(".")[0] == "matplotlib":
                self.matplotlib_first = self.matplotlib_first or node.lineno
        # get matplotlib aliases
        for child in node.names:
            if child.name == "matplotlib":
                self.matplotlib_aliases.append(child.asname or child.name)
        # get pyplot aliases
        for child in node.names:
            if child.name == "matplotlib.pyplot":
                self.pyplot_aliases.append(child.asname or child.name)
            elif child.name == "matplotlib":
                name = (child.asname or "matplotlib") + ".pyplot"
                self.pyplot_aliases.append(name)
        self.make_pyplot_paths()

    def visit_ImportFrom(self, node):
        # get line number of the first matplotlib import
        if node.module.split(".")[0] == "matplotlib":
            self.matplotlib_first = self.matplotlib_first or node.lineno
        # get pyplot aliases
        for child in node.names:
            if node.module == "matplotlib" and child.name == "pyplot":
                self.pyplot_aliases.append(child.asname or child.name)
        self.make_pyplot_paths()

    def visit_Expr(self, node):
        # get matplotlib custom options that need to be turned off
        if hasattr(node.value, "func") and hasattr(node.value.func, "value"):
            value = node.value.func.value
            # detect pyplot interactive mode
            if node.value.func.attr == "ion":
                if hasattr(value, "id") and (value.id,) in self.pyplot_paths or hasattr(
                        value.value, "id") and (value.value.id, value.attr) in self.pyplot_paths:
                    self.pyplot_interactive_linenos.append(node.lineno)
            # detect matplotlib custom backends
            if node.value.func.attr == "use":
                if hasattr(
                        value, "id") and value.id in self.matplotlib_aliases:
                    self.matplotlib_backend_linenos.append(node.lineno)
        # detect IPython magic functions: ``%matplotlib ...`` statements are
        # converted to ``get_ipython().run_line_magic('matplotlib', '...')``
        if hasattr(node.value, "func") and hasattr(node.value.func, "value") \
                and hasattr(node.value.func.value, "func") \
                and hasattr(node.value.func.value.func, "id") \
                and node.value.func.value.func.id == "get_ipython":
            self.ipython_magic_linenos.append(node.lineno)


def disable_matplotlib_gui(code):
    """
    Use the matplotlib Agg backend (no GUI).
    """
    # remove magic functions that use the percent sign syntax
    code = protect_ipython_magics(code)
    tree = ast.parse(code)
    visitor = GetMatplotlibPyplot()
    visitor.visit(tree)
    first_mpl_import = visitor.matplotlib_first
    # split lines
    lines = code.split("\n")
    mapping = delimit_statements(code)
    # list of lines to comment out
    lines_comment_out = []
    # remove magic functions
    lines_comment_out += visitor.ipython_magic_linenos
    # remove interactive mode
    lines_comment_out += visitor.pyplot_interactive_linenos
    # remove any custom backend
    lines_comment_out += visitor.matplotlib_backend_linenos
    # use the Agg backend
    if first_mpl_import:
        line = lines[first_mpl_import - 1]
        indentation = line[:len(line) - len(line.lstrip())]
        lines[first_mpl_import - 1] = indentation + \
            "import matplotlib as _mpl;_mpl.use('Agg');" + line.lstrip()
    # comment out lines
    for lineno_start in lines_comment_out:
        lineno_end = mapping[lineno_start]
        for lineno in range(lineno_start, lineno_end + 1):
            lines[lineno - 1] = "#" + lines[lineno - 1]
    code = "\n".join(lines)
    return code


class GetEspressomdVisualizerImports(ast.NodeVisitor):
    """
    Find line numbers and aliases of imported ESPResSo visualizers.
    """

    def __init__(self):
        self.visualizers = {
            "visualization",
            "visualization_opengl",
            "visualization_mayavi"}
        self.namespace_visualizers = {
            "espressomd." + x for x in self.visualizers}
        self.visu_items = {}

    def register_import(self, lineno, from_str, module_str, alias):
        if lineno not in self.visu_items:
            self.visu_items[lineno] = []
        if from_str:
            line = "from {} import {}".format(from_str, module_str)
        else:
            line = "import {}".format(module_str)
        if alias:
            line += " as {}".format(alias)
        self.visu_items[lineno].append(line)

    def visit_Import(self, node):
        # get visualizer alias
        for child in node.names:
            if child.name in self.namespace_visualizers:
                self.register_import(
                    node.lineno, None, child.name, child.asname)

    def visit_ImportFrom(self, node):
        if node.module in self.namespace_visualizers:
            for child in node.names:
                if child.name == "*":
                    raise ValueError("cannot use MagicMock() on a wildcard "
                                     "import at line {}".format(node.lineno))
                self.register_import(
                    node.lineno, node.module, child.name, child.asname)
        # get visualizer alias
        if node.module == "espressomd":
            for child in node.names:
                if child.name in self.visualizers:
                    self.register_import(
                        node.lineno, node.module, child.name, child.asname)


def mock_es_visualization(code):
    """
    Replace ``import espressomd.visualization_<backend>`` by a ``MagicMock()``
    when the visualization module is unavailable, by catching the
    ``ImportError()`` exception. Please note that ``espressomd.visualization``
    is deferring the exception, thus requiring additional checks.

    Import aliases are supported, however please don't use
    ``from espressomd.visualization import *`` because it hides the namespace
    of classes to be mocked.
    """
    # replacement template
    r_es_vis_mock = r"""
try:
    {0}{1}
except ImportError:
    from unittest.mock import MagicMock
    import espressomd
    {2} = MagicMock()
""".lstrip()

    def check_for_deferred_ImportError(line, alias):
        if "_opengl" not in line and "_mayavi" not in line:
            if "openGLLive" in line or "mayaviLive" in line:
                return """
    if hasattr({0}, 'deferred_ImportError'):
        raise {0}.deferred_ImportError""".format(alias)
            else:
                return """
    if hasattr({0}.mayaviLive, 'deferred_ImportError') or \\
       hasattr({0}.openGLLive, 'deferred_ImportError'):
        raise ImportError()""".format(alias)
        else:
            return ""

    visitor = GetEspressomdVisualizerImports()
    visitor.visit(ast.parse(protect_ipython_magics(code)))
    lines = code.split("\n")
    for lineno, imports in visitor.visu_items.items():
        line = lines[lineno - 1]
        indentation = line[:len(line) - len(line.lstrip())]
        lines[lineno - 1] = ""
        for import_str in imports:
            alias = import_str.split()[-1]
            checks = check_for_deferred_ImportError(import_str, alias)
            import_str_new = "\n".join(indentation + x for x in
                                       r_es_vis_mock.format(import_str, checks, alias).split("\n"))
            lines[lineno - 1] += import_str_new

    return "\n".join(lines)


def skip_future_imports_dependency(filepath):
    """
    If an import failed, all subsequent imports will be skipped. The
    fixture message provides the name of the module that failed.
    """
    global skip_future_imports
    if not skip_future_imports:
        module_name = os.path.splitext(os.path.basename(filepath))[0]
        assert module_name != ""
        skip_future_imports = module_name
    return unittest.skip("failed to import {}, skipping test!"
                         .format(skip_future_imports))
