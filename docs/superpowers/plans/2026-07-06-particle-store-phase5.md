# ParticleStore Migration — Phase 5 Implementation Plan (parameters eviction)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Evict all PARAMETER fields (id, mol_id, type, propagation, rotation/ext_flag bitfields, mass, rinertia, q, mu_E, dipm, gamma, gamma_rot, ext_force, ext_torque) plus the observable `dip_fld` into ParticleStore columns; move the cold PODs (swim, magnetodynamics, vs_relative) into host sidecars; `commit_particle` dies except exclusion flags (pack charge/dipm become store aliases); the PROPRTS ghost group goes value-based.

**Architecture:** Fifth instantiation of the proven template (audit → columns+carriers → ghost value path → flip → checkpoint). New elements: scalar columns of mixed types (int, uint8_t, double; scalar-or-Vector3d gamma per ESPRESSO_PARTICLE_ANISOTROPY via a typedef), constexpr-when-disabled accessor preservation (mass/q/rinertia/rotation/ext_flag keep their static fallbacks when the feature is off — no column exists then), host sidecar vectors for the three cold PODs (indexed by store row, rebuilt with the store, values ferried through PROPRTS + migration carriers), and the type/charge mutation surface of reaction methods (the primary physics gate).

**Spec:** section 3 phase 5; budget RE-TIGHTENS to ≤5% on ALL configs at this checkpoint (amendment expires). **Exploration report with per-field verdicts: this plan's basis (see `.superpowers/sdd/` phase-5 exploration).**

## Global Constraints

- All standing constraints (worktree, `make -j8`, identity gate, clean-core-rebuild after ABI changes — EIGHT incidents to date, controller-only benchmarks, noise allowlist, commit trailer).
- NEW standing gate since phase 4: CI-mirror (`build-ci-mirror/`, Debug + ESPRESSO_WARNINGS_ARE_ERRORS=ON) core compile green at every task that touches headers; GitHub Actions must be green after the phase push (`gh run list --repo RudolfWeeber/espresso --branch worktree-eliminate_Particle`).
- Physics gates MUST include reaction methods (`ctest -R "reaction|constant_pH|widom"` — type/charge mutations), plus the standing lees_edwards/hybrid/virtual/dpd/lb set.
- id() is written only by particle_node (creation/id-reassignment); type/q are written by ReactionAlgorithm mid-run — `on_particle_type_change` bookkeeping fires at the SETTER layer (particle_node), unaffected by accessor redirection: VERIFY, don't assume.
- Constexpr-when-disabled fields: accessors keep the static path under `#else`; columns/carriers/ghost fields are ifdef-gated identically to today's members.
- `m_particle_index` (Particle*) SURVIVES until phase 7; the id column serves kernels/ghost packing; `get_local_particle` unchanged.

## Tasks

### Task 1: Preparation audit
Grep-driven, zero behavior change: template-arithmetic/reference-binding uses of ext_force(), ext_torque(), rinertia(), mu_E(), gamma()/gamma_rot() (Vector3d-typed per config!), dip_fld() (accumulation sites in magnetostatics + GpuParticleData `dip_fld()[j] +=`), and any `auto &x = p.<param>()` bindings (scalar refs like `auto &t = p.type()` — convert to flip-surviving forms). Scalar reads (mass/q/dipm/type/id in arithmetic) are proxy-free post-flip (element references) — verify the planned accessor shape: scalar columns return real `T&` (like lees_edwards_offset), so MOST scalar sites need nothing; Vector3d params get VectorReference treatment. vs_relative/swim/magnetodynamics accessor uses: inventory for the sidecar API (return references to sidecar structs — semantics unchanged). Gates: unit, physics battery incl. reaction tests, identity. Commit.

### Task 2: Columns, sidecars, carriers
- Scalar columns: id/mol_id/type/propagation (`DualView<int*>`), rotation/ext_flag (`DualView<uint8_t*>` — ifdef), mass/q/dipm (`DualView<double*>` — ifdef). Vector columns: rinertia/mu_E/ext_force/ext_torque (ifdef), dip_fld (observable category — zeroed via init like force). Gamma: `using GammaValue = double | Utils::Vector3d` per ESPRESSO_PARTICLE_ANISOTROPY; column typedef follows.
- Host sidecars: `std::vector<ParticleParametersSwimming>`, `std::vector<ThermalStonerWohlfarthParameters>`, `std::vector<VirtualSitesRelativeParameters>` (ifdef-gated), sized/rebuilt with the store (preserve-by-old-row in assign_row, seed from carriers).
- DORMANT carriers for everything (the three PODs ride as whole-POD carriers).
- Defaults must match current member defaults EXACTLY (id -1, mol_id 0, type 0, propagation SYSTEM_DEFAULT, bitfields 0, mass 1.0, rinertia {1,1,1}, q 0, dipm 0, gamma -1 or {-1,-1,-1}, vs_relative quats identity...).
- Unit tests per established pattern. ABI → clean rebuild + CI-mirror compile. Commit.

### Task 3: Ghost PROPRTS value path
serialize_and_reduce PROPERTIES branch → direction-first value-based, preserving the EXACT field order/ifdef structure (wire byte-identical; ~89-220 B/particle); ParameterRowContext for the hot scalar/vector columns + sidecar lookups for PODs; calc_transmit_size masking (PROPRTS is the last Particle{}-sized group — after this, calc_transmit_size may become fully constant: check what remains (BONDS/RATTLE) and simplify if trivial). PROPRTS is resort-only — perf-light, correctness-heavy. Four-combination enumeration mandatory. Gates + identity. Commit.

### Task 4: THE FLIP
- ParticleProperties member `p` leaves Particle (struct fully hollowed except BondList/exclusions/local flags — check what remains: bl, el, rattle, ghost flag...); accessors: scalar params → `T&`/value via column element refs with constexpr-when-disabled preserved and carrier fallback when detached; Vector3d params → VectorReference/value; PODs → sidecar references (attached) / carrier (detached).
- Serialize: `ar & p` leaves; parameter carriers added (once).
- Pack: charge/dipm views become store-column ALIASES (bind_pack_store_views pattern); id/type/mass pack views likewise aliased or kernels read columns via row translation — pick the simpler consistent form; commit_particle reduces to exclusion flags only.
- Python: setters route via accessors (verify type/q setter → on_particle_type_change path & reaction methods); getters via snapshot store (verify remote reads).
- Sync-point audit for parameter readers post-topology-change (analysis by type, reaction loops).
- Fixtures sweep; serialization test expectations.
- Gates: full battery + REACTION tests + identity + CI-mirror. Commit.

### Task 5: Checkpoint + ship
Guard; unit; hybrid until-fail:3; CI-mirror; full check_python(+skip_long); identity; HDF5 side-build; CONTROLLER benchmark gate 3 reps at **restored ≤5% ALL configs** (charge/dipm copy deletion + commit_particle death should help lj1; if any config fails 5%, HALT with the table — the amendment has expired). Spec adjudication updates; push; **verify GitHub CI green**; ledger+memory.
