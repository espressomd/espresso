#!/usr/bin/env bash
# Guard for the ParticleStore migration (see docs/superpowers/specs/
# 2026-07-03-array-based-particle-storage-design.md, section 3):
# cell particle storage may only be mutated through the primitives in
# src/core/cell_system/ParticleListOperations.hpp. This script fails if a
# direct mutation pattern reappears in src/core (unit tests excluded:
# they construct fixtures directly).
#
# Phase 7c: a Cell's committed content collapsed from a Bag<int> row bag to a
# contiguous (offset, count) store-row RANGE (Cell::set_range) plus a staging
# buffer of not-yet-committed row references (Cell::staged()). A row dropped
# mid-epoch is marked pending-removed on the store
# (ParticleStore::mark_pending_removal). The choke points now: mark a range
# position pending-removed on drop, collapse the range to empty on clear/ghost
# resize, and stage inserts/ghost slots. The store rebuild writes the range back
# (set_range). The out-of-band-mutation hazards this guard catches are therefore
# a raw set_range / mark_pending_removal / staged()-mutation OUTSIDE the choke
# points + the rebuild (the range collapse choke points replace the pre-7c
# rows()-bag mutations).
#
# The scan is delegated to an embedded Python program because the mutation
# call can be split across several physical lines by clang-format, e.g.
#     object.get_local_cells()[index]
#         ->rows()
#         .insert(row);
# A line-based grep misses this; the Python scanner reads each file as one
# string and matches the pattern across newlines.
#
# ALIAS LIMITATION: this check cannot see mutations performed through a
# local reference/alias to a cell's storage or a hoisted store handle, e.g.
#     auto &st = cell->store();
#     st.mark_pending_removal(row);
# because the mutating call no longer names the guarded surface. A green run is
# a tripwire that catches the common, direct form; it is NOT a proof that no
# direct cell-storage mutation exists. Reviewers must still watch for aliased
# mutations by hand.
#
# STORE-DIRECT WRITES OUT OF SCOPE BY DESIGN: this guard protects the
# cell-storage choke points (rows()/staged()/particles() mutations).  Direct
# writes to the ParticleStore itself -- e.g. store.assign_row(...) or column
# writes through a view -- are NOT checked here.  In practice assign_row is
# only callable inside a begin_rebuild/finish_rebuild bracket and column writes
# go through Particle views that are only valid during a live generation, so
# the risk is low; but the "single choke point" guarantee extends to cell
# storage only, not to the underlying store columns.
#
# EXCEPTIONS: (1) ParticleListOperations.hpp defines the choke points (drop_row
# marks pending-removal; clear/ghost-resize collapse the range). (2)
# CellStructure.cpp's ensure_particle_store_synchronized owns the store rebuild
# and legitimately writes every cell's range back (set_range) as it renumbers
# rows. (3) Cell.hpp DEFINES set_range (the setter itself, not a mutation site).
# All three files are excluded below.
set -u
cd "$(dirname "$0")/../.."
exec python3 - "$@" << 'PYEOF'
import os
import re
import sys

ROOT = "src/core"
EXCLUDED_DIRECTORIES = ("src/core/unit_tests/",)
# Choke points (ParticleListOperations.hpp) legitimately mark rows
# pending-removed and collapse cell ranges; the store rebuild in
# CellStructure.cpp legitimately writes every cell's range back (set_range) as
# it renumbers rows; Cell.hpp defines set_range. All are the sanctioned owners
# of cell storage and are exempt.
EXCLUDED_FILES = (
    "src/core/cell_system/ParticleListOperations.hpp",
    "src/core/cell_system/CellStructure.cpp",
    "src/core/cell_system/Cell.hpp",
)

# Phase 7c: a Cell's committed content is a contiguous (offset, count) store-row
# range (written back by set_range) plus a staging buffer (Cell::staged()); a
# mid-epoch drop marks the store row pending-removed
# (ParticleStore::mark_pending_removal). Direct use of any of these -- outside
# the CellParticleStorage choke points and the store rebuild -- is the
# out-of-band-mutation hazard. Match:
#   - <cell>.set_range( ... )              (raw range writeback)
#   - <store>.mark_pending_removal( ... )  (raw drop outside drop_row)
#   - staged()->/.<mutator>( ... )         (raw staging mutation)
# allowing arbitrary whitespace (including newlines) around the '.'/'->'
# operator so clang-format-wrapped calls are caught.
PATTERN = re.compile(
    r"(?:"
    r"(?:->|\.)\s*(?:set_range|mark_pending_removal)\s*\("
    r"|"
    r"staged\(\)\s*(?:->|\.)\s*"
    r"(?:insert|erase|clear|resize|push_back|emplace_back|pop_back|emplace)"
    r"\s*\("
    r")")


def normalized(path):
    return path.replace(os.sep, "/")


def offending_files():
    findings = []
    for dirpath, _, filenames in os.walk(ROOT):
        for filename in filenames:
            if not filename.endswith((".cpp", ".hpp")):
                continue
            path = normalized(os.path.join(dirpath, filename))
            if path.startswith(EXCLUDED_DIRECTORIES):
                continue
            if path in EXCLUDED_FILES:
                continue
            with open(path, encoding="utf-8", errors="replace") as handle:
                text = handle.read()
            for match in PATTERN.finditer(text):
                line_number = text.count("\n", 0, match.start()) + 1
                snippet = " ".join(match.group(0).split())
                findings.append((path, line_number, snippet))
    return findings


def main():
    findings = offending_files()
    if findings:
        print("ERROR: direct cell-storage mutation outside "
              "CellParticleStorage:")
        for path, line_number, snippet in findings:
            print(f"  {path}:{line_number}: {snippet}")
        return 1
    print("OK: no direct cell-storage mutations found.")
    return 0


sys.exit(main())
PYEOF
