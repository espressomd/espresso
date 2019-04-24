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
import sys
import unittest
import importlib
import espressomd
if sys.version_info >= (3, 3):
    from unittest.mock import MagicMock
else:
    from mock import MagicMock


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
    - use random seeds for the RNG in NumPy and ESPResSo
    - temporarily move to the directory where the script is located

    :param filepath: python script to import
    :type  filepath: str
    :param gpu: whether GPU is necessary or not
    :type  gpu: bool
    :param substitutions: custom text replacement operation (useful to edit out
       calls to the OpenGL or Mayavi visualizers' :meth:`run` method)
    :type  substitutions: function
    :param cmd_arguments: command line arguments, i.e. sys.argv without the
       script path
    :type  cmd_arguments: list
    :param script_suffix: suffix to append to the configured script (useful
       when a single module is being tested by multiple tests in parallel)
    :type  script_suffix: str
    :param random_seeds: if ``True``, use random seeds in RNGs
    :type  random_seeds: bool
    :param mock_visualizers: if ``True``, substitute ES visualizers with
       `Mock()` classes in case of `ImportError()` (set to ``False`` if an
       `ImportError()` is relevant to your test)
    :type  mock_visualizers: bool
    :param move_to_script_dir: if ``True``, move to the script's directory
       (useful when the script needs to load files hardcoded as relative paths,
       or when files are generated and need cleanup); this is enabled by default
    :type  move_to_script_dir: bool
    :param \*\*parameters: global variables to replace
    :type  \*\*parameters: int, float, bool
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
    # use random seeds for ES and NumPy RNGs
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


def set_cmd(code, filepath, cmd_arguments):
    assert isinstance(cmd_arguments, list) \
        or isinstance(cmd_arguments, tuple)
    sys_argv = list(map(str, cmd_arguments))
    sys_argv.insert(0, os.path.basename(filepath))
    re_import_sys = re.compile("^import[\t\ ]+sys[\t\ ]*$", re.M)
    re_import_argparse = re.compile("^import[\t\ ]+argparse[\t\ ]*$", re.M)
    if re_import_sys.search(code) is not None:
        code = re_import_sys.sub("\g<0>\nsys.argv = " + str(sys_argv), code, 1)
    elif re_import_argparse.search(code) is not None:
        code = re_import_argparse.sub("\g<0>\nimport sys\nsys.argv = "
                                      + str(sys_argv), code, 1)
    else:
        raise AssertionError("module sys (or argparse) is not imported")
    old_sys_argv = list(sys.argv)
    return code, old_sys_argv


def substitute_variable_values(code, **parameters):
    for variable, value in parameters.items():
        assert variable in code, "variable {} not found".format(variable)
        re_var = re.compile("^(\t|\ {,4})(" + variable + ")(?= *=[^=])", re.M)
        assert re_var.search(code) is not None, \
            "variable {} has no assignment".format(variable)
        code = re_var.sub(
            r"\g<1>\g<2> = " + repr(value) + r"; _\g<2>__original", code)
    return code


def set_random_seeds(code):
    # delete explicit ESPResSo seed
    aliases = re.findall(r"([^\s;]+) *= *(?:espressomd\.)?System *\(", code)
    pattern = r"(?<=[\s;]){}\.(?:seed|random_number_generator_state)(?= *=[^=])"
    subst = "{}.set_random_state_PRNG(); _random_seed_es__original"
    for varname in set(aliases):
        code = re.sub(pattern.format(varname), subst.format(varname), code)
    # delete explicit NumPy seed
    code = re.sub(r"(?<=[\s;])(?:numpy|np)\.random\.seed *(?=\()",
                  "_random_seed_np = (lambda *args, **kwargs: None)", code)
    return code


def disable_matplotlib_gui(code):
    """
    Use the matplotlib Agg backend (no GUI).
    """
    # find under which name matplotlib was imported
    re_mpl_aliases = [
        re.compile(r"^[\t\ ]*import[\t\ ]+(matplotlib)[\t\ ]*$", re.M),
        re.compile(r"^[\t\ ]*import[\t\ ]+matplotlib[\t\ ]+as[\t\ ]+([^\s;]+)",
                   re.M)]
    aliases = set(x for re_mpl in re_mpl_aliases for x in re_mpl.findall(code))
    # remove any custom backend
    for alias in aliases:
        code = re.sub(r"^[\t\ ]*" + alias + r"\.use\(([\"']+).+?\1[\t\ ]*\)",
                      "", code, 0, re.M)
    # use the Agg backend
    code = re.sub(r"^([\t\ ]*)(?=(?:from|import)[\t\ ]+matplotlib[\.\s])",
                  r"\g<1>import matplotlib as _mpl;_mpl.use('Agg');",
                  code, 1, re.M)
    return code


def mock_es_visualization(code):
    """
    Replace `import espressomd.visualization_<backend>` by a `MagicMock()` when
    the visualization module is not installed, by catching the `ImportError()`
    exception. Please note that `espressomd.visualization` is deferring the
    exception, thus requiring additional checks. Import aliases are supported,
    however please don't use `from espressomd.visualization import *` because
    it hides the namespace of classes to be mocked.
    """
    # consider all legal import statements in Python3
    # (the ordering follows regex precedence rules)
    re_es_vis_import = re.compile(r"""
 ^from\ espressomd\ import\ (?:visualization(?:_opengl|_mayavi)?)\ as\ (\S+)
|^from\ espressomd\ import\ (visualization(?:_opengl|_mayavi)?)
|^from\ espressomd\.visualization(?:_opengl|_mayavi)?\ import\ \S+\ as\ (\S+)
|^from\ espressomd\.visualization(?:_opengl|_mayavi)?\ import\ (\S+)
|^import\ espressomd\.visualization(?:_opengl|_mayavi)?\ as\ (\S+)
|^import\ (espressomd\.visualization(?:_opengl|_mayavi)?)
""".replace(r"\ ", r"[\t\ ]+"), re.VERBOSE | re.M)
    # replacement template
    r_es_vis_mock = r"""
try:
    {0}{1}
except ImportError:
    from {2} import MagicMock
    import espressomd
    {3} = MagicMock()
""".lstrip()
    # cannot handle "from espressomd.visualization import *"
    re_es_vis_import_namespace = re.compile(
        r"^from\ espressomd\.visualization(?:_opengl|_mayavi)?\ import\ \*"
        .replace(r"\ ", r"[\t\ ]+"), re.M)
    m = re_es_vis_import_namespace.search(code)
    assert m is None, "cannot use MagicMock() at line '" + m.group(0) + "'"

    def check_for_deferred_ImportError(s, alias):
        if "_opengl" not in s and "_mayavi" not in s:
            if "openGLLive" in s or "mayaviLive" in s:
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

    def substitution_es_vis_import(m):
        mock_module = "unittest.mock" if sys.version_info >= (3, 3) else "mock"
        alias = [x for x in m.groups() if x is not None][0]
        checks = check_for_deferred_ImportError(m.group(0), alias)
        return r_es_vis_mock.format(m.group(0), checks, mock_module, alias)

    # handle deferred ImportError
    code = re_es_vis_import.sub(substitution_es_vis_import, code)
    return code


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
