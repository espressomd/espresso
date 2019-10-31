import os
import ast
import textwrap

samples_dir = '@CMAKE_SOURCE_DIR@/samples/'


def get_docstring(filenames):
    docstrings = {}
    for filename in filenames:
        with open(samples_dir + filename) as f:
            src = f.read()
        tree = ast.parse(src)
        module = next(ast.walk(tree))
        docstrings[filename] = ast.get_docstring(module)
    return docstrings


# extract docstrings
samples = [x for x in os.listdir(samples_dir) if x.endswith('.py')]
samples += ['immersed_boundary/sampleImmersedBoundary.py',
            'object_in_fluid/motivation.py']
docstrings = get_docstring(samples)


# write documentation
sphinx_tpl = '* :file:`{}`\n{}\n'
with open('@CMAKE_CURRENT_BINARY_DIR@/samples.rst', 'w') as f:
    for filename in sorted(docstrings, key=lambda x: x.lower()):
        docstring = (docstrings[filename] or '').replace('ESPResSo', '|es|')
        paragraph = textwrap.indent(docstring, prefix='    ')
        f.write(sphinx_tpl.format(filename, paragraph))
