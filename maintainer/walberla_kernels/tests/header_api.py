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
Extraction of the function signatures of a C++ header that defines
namespace-level ``inline`` functions or declares namespace-level functions
(``declarations=True``), used to compare the generated field
accessors with a reference header. Type spellings are canonicalized, such that
signatures that differ only in the namespace qualifiers or the spelling of
integer literals (``uint_t{19u}``, ``19u``, ``19``) are equal.
"""

import collections
import json
import re


def strip_comments_and_strings(text):
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
    text = re.sub(r"//[^\n]*", "", text)
    return re.sub(r'"(?:\\.|[^"\\])*"', '""', text)


def split_top_level(text, sep=","):
    """Split at the separators that are not enclosed in brackets."""
    parts, depth, current = [], 0, ""
    for ch in text:
        if ch in "<([{":
            depth += 1
        elif ch in ">)]}":
            depth -= 1
        if ch == sep and depth == 0:
            parts.append(current)
            current = ""
        else:
            current += ch
    if current.strip():
        parts.append(current)
    return [p.strip() for p in parts]


def canonical_type(type_):
    t = " ".join(type_.split())
    t = t.replace("walberla::field::", "").replace("walberla::", "")
    t = re.sub(r"\bfield::", "", t)
    t = re.sub(r"uint_t\{\s*(\d+)u\s*\}", r"\1", t)
    t = re.sub(r"\b(\d+)u\b", r"\1", t)
    t = re.sub(r"\s*([<>,&*])\s*", r"\1", t)
    t = t.replace(">>", "> >")
    base, declarator = re.match(r"^(.*?)([&*]*)(?:\s*const)?$", t).groups()
    base = base.strip()
    is_const = False
    if re.match(r"^const\b", base):
        is_const = True
        base = re.sub(r"^const\s*", "", base)
    if re.search(r"\bconst$", base):
        is_const = True
        base = re.sub(r"\s*\bconst$", "", base)
    if not declarator:
        # top-level const of a by-value parameter is not part of the API
        is_const = False
    return f"{base}{' const' if is_const else ''}{declarator}"


def canonical_expression(expr):
    return re.sub(r"\s+", "", expr.replace("walberla::", ""))


def parse_parameter(text):
    """Return the tuple (type, name, default argument) of a parameter."""
    default = None
    if "=" in text:
        text, default = text.split("=", 1)
        default = canonical_expression(default)
    type_, name = re.match(
        r"^(.*?)(\w+)\s*$", text.strip(), flags=re.S).groups()
    return canonical_type(type_), name, default


def _declaration(statement):
    """
    Return ``(return type, name, parameter text)`` of a function declaration
    ``<return type> <name>(<parameters>)``, or None for any other statement.
    """
    match = re.match(r"^(.*?)\b(\w+)\s*\((.*)\)\s*$", statement, flags=re.S)
    if match is None or not match.group(1).strip():
        return None
    if re.match(r"^\s*(using|typedef|static_assert)\b", statement):
        return None
    return match.group(1), match.group(2), match.group(3)


def functions(text, declarations=False, include_kernels=False):
    """
    Return the tuples (namespace, name, parameters, return type) of all
    namespace-level ``inline`` functions of a header, where the parameters
    are tuples (type, name, default argument). With ``declarations=True``,
    the namespace-level function declarations (no definitions) are returned
    instead. The functions of the namespaces ``kernels`` (the pystencils
    kernels of a header-only file) are not part of the API, unless
    ``include_kernels`` is set.
    """
    text = strip_comments_and_strings(text)
    if declarations:
        # preprocessor directives, including line continuations
        text = re.sub(r"^[ \t]*#(?:[^\n]|\\\n)*$", "", text, flags=re.M)
    signature = re.compile(r"inline\s+(.*?)\b(\w+)\s*\(", re.S)
    stack = []  # namespace names, or None for any other block
    result = []
    last_boundary = 0
    i = 0

    def namespace_path():
        return "::".join(stack).replace("walberla::", "", 1)

    def is_api():
        return include_kernels or stack[-1:] != ["kernels"]

    while i < len(text):
        ch = text[i]
        if ch == "{":
            match = re.search(r"namespace\s+([\w:]+)\s*$",
                              text[last_boundary:i])
            stack.append(match.group(1) if match else None)
            last_boundary = i + 1
        elif ch == "}":
            if stack:
                stack.pop()
            last_boundary = i + 1
        elif ch == ";":
            if (declarations and all(s is not None for s in stack)
                    and is_api()):
                declaration = _declaration(text[last_boundary:i])
                if declaration is not None:
                    return_type, name, params = declaration
                    result.append((
                        namespace_path(), name,
                        tuple(parse_parameter(p)
                              for p in split_top_level(params)),
                        canonical_type(return_type)))
            last_boundary = i + 1
        elif (not declarations and text.startswith("inline", i)
              and (i == 0 or not text[i - 1].isalnum())
              and all(s is not None for s in stack)):
            match = signature.match(text, i)
            if match:
                depth, end = 1, match.end()
                while depth:
                    depth += {"(": 1, ")": -1}.get(text[end], 0)
                    end += 1
                params = tuple(
                    parse_parameter(p)
                    for p in split_top_level(text[match.end():end - 1]))
                if is_api():
                    result.append((namespace_path(), match.group(2), params,
                                   canonical_type(match.group(1))))
                i = end
                continue
        i += 1
    return result


# the keys under which two lists of functions are compared, from the coarsest
# to the finest
SIGNATURE_KEYS = {
    "parameter types": lambda f: (f[0], f[1], tuple(p[0] for p in f[2])),
    "parameter names": lambda f: (f[0], f[1],
                                  tuple((p[0], p[1]) for p in f[2])),
    "default arguments": lambda f: (f[0], f[1], f[2]),
    "return type": lambda f: f,
}


def signature_counters(text, declarations=False):
    """Multisets of the functions of a header, for each key."""
    return counters(functions(text, declarations))


def to_records(funcs):
    """Sorted, JSON-serializable list of records of a list of functions."""
    records = [{"namespace": f[0], "name": f[1],
                "parameter types": [p[0] for p in f[2]],
                "parameter names": [p[1] for p in f[2]],
                "defaults": [p[2] for p in f[2]],
                "return type": f[3]} for f in funcs]
    return sorted(records, key=lambda r: json.dumps(r, sort_keys=True))


def from_records(records):
    """Inverse of `to_records`."""
    return [(r["namespace"], r["name"],
             tuple(zip(r["parameter types"], r["parameter names"],
                       r["defaults"])), r["return type"])
            for r in records]


def load_records(path):
    with open(path) as f:
        return from_records(json.load(f))


def counters(funcs):
    """Multisets of a list of functions, for each key."""
    return {label: collections.Counter(key(f) for f in funcs)
            for label, key in SIGNATURE_KEYS.items()}
