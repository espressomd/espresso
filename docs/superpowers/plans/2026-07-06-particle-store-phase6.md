# ParticleStore Migration — Phase 6 Implementation Plan (ragged data eviction)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Evict the last owned non-POD members of `Particle` — `BondList bl` and exclusions `el` — into ParticleStore host sidecars indexed by store row, and evict `ParticleRattle::correction` as a plain observable column; after this phase `Particle` holds only the store attach handle, the migration carriers (die in phase 7), and `ParticleLocal`.

**Architecture:** Sixth instantiation of the proven template, but LEANER: the exploration report (`.superpowers/sdd/phase6-exploration.md` — this plan's basis, read it before each task) established that every bond/exclusion consumer (LocalBondState flattening in `CellStructure.cpp:237-360`, the three bonded Kokkos kernels, the `has_exclusion` pack write, `do_nonbonded`, the ghost BONDS buffer) reads through `p.bonds()`/`p.exclusions()`, so the flip is a single-point accessor redirect. Ragged carriers use the DUAL-ROLE pattern: the struct members stay as the detached/migration storage, the attached accessor redirects to the sidecar. The ghost BONDS wire path (separate `CommBuf::bondbuf`, boost binary archive, excluded from `calc_transmit_size`) is kept EXACTLY as is — do NOT columnarize bonds.

**Spec:** section 3 phase 6. Perf: phase 6 is perf-neutral on the gated benchmarks by construction (no bonds/exclusions in those configs); the checkpoint still runs the full gate at the AMENDED budget (2026-07-06 amendment: 1000ppc configs ≤12%, 4000ppc configs ≤5%).

## Global Constraints

- All standing constraints: work in this worktree only; `make -j8` (never `-j$(nproc)`); identity gate both modes vs `.superpowers/sdd/phase1-identity-reference.txt`; controller-only benchmarks; commit trailer `Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>`.
- Stale-ABI rule (NINE incidents): after ANY change to `Particle.hpp`/`ParticleStore.hpp`/pack headers, clean-rebuild `espresso_core` AND `make -C build -j8 unit_tests_executables` before trusting any test result.
- CI gates at every header-touching task: CI-mirror core compile (`build-ci-mirror/` — symlink to `/ssd/weeber/eliminate_Particle_builds/`, Debug+ESPRESSO_WARNINGS_ARE_ERRORS=ON) 0 errors/0 warnings; maxset core compile (`build-maxset/`, all ifdef features) 0 errors. GitHub Actions must be green after the phase push (`gh run list --repo RudolfWeeber/espresso --branch worktree-eliminate_Particle`).
- Physics gates MUST include the bond surface: `ctest -R "collision|bond|rigid|exclusion|virtual|lees_edwards|hybrid|dpd"` plus the standing set. `rigid_bond` is the primary RATTLE gate; `collision_detection` tests are the primary mid-run bond-mutation gate.
- Test patterns for any ghost/storage change include `lees_edwards`, `hybrid_decomposition`, `virtual_sites`.
- Full-word names over abbreviations in new identifiers.
- The ghost BONDS wire format and `calc_transmit_size` bond exclusion are behavior-frozen this phase.

## Tasks

### Task 1: Preparation audit
Grep-driven, zero behavior change. Inventory every `p.bonds()`/`p.exclusions()`/`p.rattle_params()` (or equivalent rattle accessor) site and classify: (a) reference bindings `auto &b = p.bonds()` — these SURVIVE the flip (sidecar accessors return real `BondList&`, not proxies) but verify none outlive a store rebuild (dangling across `assign_row`/resort); (b) mutation sites on attached particles mid-run (collision detection, bond breakage, reaction deletion) — confirm each is followed by the existing resort/dirty machinery, no new hooks expected; (c) `Particle::serialize` bond/exclusion legs — confirm the envelope form. Rattle: inventory `correction` reads/writes (`rattle.cpp`, ghost reduction `ghosts.cpp:1838-1852`). Fix only what would break under the flip (per exploration: expected near-zero). Gates: unit, bond-focused physics battery, identity. Commit `core: prepare ragged-data call sites for sidecar accessors (phase 6)`.

### Task 2: Sidecars, column, carriers (dormant)
- ParticleStore ragged host sidecars: `std::vector<BondList> m_bonds` and (ESPRESSO_EXCLUSIONS) `std::vector<ExclusionList> m_exclusions` (use the exact current type of `Particle::el`), following the phase-5 POD-sidecar machinery verbatim: sized in `begin_rebuild` (old vector swapped), `preserve_or_seed_sidecar` in `assign_row` (move-from-old-row where possible — ragged sizes make copies expensive; use `std::move` on the old element, matching vector element move semantics), cleared in `release_columns`; accessors `bonds_sidecar_reference(row)` / `exclusions_sidecar_reference(row)` via the existing `sidecar_reference` helper. Ghost rows get default-empty elements (exploration finding 4: ghost exclusion sidecars stay empty; ghost BOND sidecars are filled only by the ghost BONDS unpack path post-flip).
- RATTLE observable column: `Column m_rattle_correction` (Vector3d, both generations, observable category like force — no carrier, no seed-from-carrier; preserve-by-row like force for mid-iteration resorts is NOT needed since SHAKE zeroes it each iteration, but keep preserve-or-default for uniformity with dip_fld).
- DORMANT: nothing reads the sidecars/column yet; Particle unchanged.
- Unit tests per the established preservation/seeding pattern: bonds sidecar preserve across rebuild (non-empty BondList moves intact), exclusions sidecar, rattle column preserve/default. ABI → clean rebuild incl. unit_tests_executables + CI-mirror + maxset compiles. Commit `core: ParticleStore ragged sidecars and rattle column (phase 6)`.

### Task 3: THE FLIP
- `Particle`: `bl`/`el` members RENAMED to `m_migration_bonds`/`m_migration_exclusions` (dual-role: detached storage AND migration envelope — exploration finding 3; do NOT tag bitwise-serializable). Accessors `bonds()`/`exclusions()` (both constnesses): attached → sidecar element reference, detached → member. Detach path (`~line 1281` pattern) syncs `m_migration_bonds = detached_bonds()` etc. alongside the existing carriers. `ParticleRattle::correction` → column accessor (attached-only assert like `force()`); `rattle` struct member shrinks or dies accordingly (its remaining field(s) — check exploration §6 — stay as-is if any).
- `Particle::serialize`: bond/exclusion legs now serialize the migration members (same boost form → wire/fetch-cache/checkpoint semantics unchanged; the remote `get_particle_data` snapshot path must still deliver bonds — verify per exploration finding 8).
- Ghost BONDS branch: NO wire change — it packs via `p.bonds()` which now reads the sidecar for attached ghosts; verify the unpack side writes through the accessor into ghost sidecar rows, and the rebuild-before-bond-unpack ordering (`update_ghosts_and_resort_particle`) still holds so ghost rows exist when the unpack runs.
- `commit_particle` exclusion-flag write, LocalBondState flattening, bonded kernels: untouched by construction — verify by grep that none bypass the accessors.
- Sync-point check: bond mutations on attached particles (collision detection binds mid-run) write the sidecar through the accessor; confirm the existing resort/Verlet-rebuild triggers still fire (no pack coherence issue — bonds are not packed).
- Fixtures sweep: detached `Particle{}` construction + bond add in tests; serialization test expectations (size unchanged — form is identical).
- Gates: unit ALL; physics battery incl. collision/bond/rigid_bond/exclusion + standing set; identity BOTH modes; CI-mirror + maxset. Commit `core: evict bonds, exclusions and rattle correction into ParticleStore (phase 6)`.

### Task 4: Checkpoint + ship
Guard script; unit; hybrid until-fail:3; CI-mirror FULL build; full check_python(+skip_long); identity; HDF5 side-build; CONTROLLER benchmark gate 3 reps at the AMENDED budget (1000ppc ≤12%, 4000ppc ≤5%) — expected neutral; investigate anything beyond noise vs the phase-5 gate CSV (`.superpowers/sdd/phase5-final-gate.csv`). Spec adjudication updates; push to `me`; verify GitHub CI green; ledger+memory.
