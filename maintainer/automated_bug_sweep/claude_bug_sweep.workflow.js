// Claude Code Workflow script (NOT a standalone node script).
// Run it from a Claude Code session anywhere inside the checkout with:
//   Workflow({ scriptPath: "maintainer/core_bug_sweep.workflow.js" })
// Optional args:
//   { root: "/abs/path/to/repo" }      -> skip git-root discovery, use this checkout
//   { only: ["electrostatics","cell_system"] } -> restrict to units whose name contains
//                                                 any of these tokens (pilot / re-run a slice)
//   { singleVerifier: true }           -> 1 verifier per finding instead of 2 for high/critical
//   { maxLinesPerUnit: 3500 }          -> tune how aggressively big directories are split
//   { includeGenerated: true }         -> also scan machine-generated source (.../generated_kernels);
//                                         skipped by default as low-signal + large agent cost
//
// The workflow is REUSABLE and self-scoping:
//   * It discovers the repository ROOT itself (git rev-parse), so nothing is hard-coded.
//   * It enumerates ALL source under src/ (core + utils + shapes + particle_observables +
//     script_interface + python + walberla_bridge/scafacos glue + ...), EXCLUDING tests,
//     then partitions the file list into review units deterministically in JS.
//   * It builds python -> script_interface -> core CALL-PATH units to audit the public-API seam.
//   * Every finding is adversarially re-verified (refute-by-default) and annotated with a
//     REACHABILITY analysis (does any calling code exclude the triggering case?).
// External library internals (walberla, heffte, kokkos/cabana, boost, MPI, ScaFaCoS) are NOT
// reviewed; reviewers inspect only ESPResSo's own logic and how it calls those libraries.
// Returns { stats, survivors, synth }.

export const meta = {
  name: 'src-bug-sweep',
  description: 'Self-scoping adversarial bug sweep of all of ESPResSo src/ (+ python->script_interface->core call paths)',
  phases: [
    { title: 'Scope', detail: 'discover repo root + enumerate src/ source files (no tests, no vendored libs)' },
    { title: 'Find', detail: 'one bug-hunter per review unit + python->SI->core call-path reviewers' },
    { title: 'Verify', detail: 'adversarial refute-by-default verifiers with reachability analysis, pipelined per finding' },
    { title: 'Synthesize', detail: 'dedupe + prioritize (silently-wrong-physics first; latent/unreachable last)' },
  ],
}

// ---- schemas -------------------------------------------------------------
const SEVERITY = ['critical', 'high', 'medium', 'low']
const CATEGORY = ['memory_safety', 'integer_overflow', 'concurrency_mpi', 'numerics', 'logic', 'resource_leak', 'api_contract', 'md_invariant', 'call_path', 'other']
const REACH_STATUS = ['unreachable', 'partially_reachable', 'reachable', 'unknown']

const SCOPE_SCHEMA = {
  type: 'object', additionalProperties: false,
  required: ['root', 'files', 'file_count', 'total_lines', 'excluded_note'],
  properties: {
    root: { type: 'string', description: 'absolute repository root (output of git rev-parse --show-toplevel)' },
    files: {
      type: 'array', description: 'every C/C++ source file under src/ (excluding tests), repo-relative',
      items: {
        type: 'object', additionalProperties: false, required: ['path', 'lines'],
        properties: {
          path: { type: 'string', description: 'repo-relative path, e.g. src/core/electrostatics/p3m.cpp' },
          lines: { type: 'integer', description: 'line count (wc -l)' },
        },
      },
    },
    python_modules: {
      type: 'array', description: 'repo-relative paths of src/python espressomd modules (*.py/*.pyx) for call-path review',
      items: { type: 'string' },
    },
    file_count: { type: 'integer', description: 'length of files[] (self-check)' },
    total_lines: { type: 'integer', description: 'sum of lines over files[] (self-check)' },
    excluded_note: { type: 'string', description: 'what was excluded (test dirs, etc.) and roughly how many' },
  },
}

const FIND_SCHEMA = {
  type: 'object', additionalProperties: false, required: ['findings'],
  properties: { findings: { type: 'array', items: {
    type: 'object', additionalProperties: false,
    required: ['title', 'file', 'severity', 'category', 'trigger_path', 'impact', 'evidence', 'confidence'],
    properties: {
      title: { type: 'string' },
      file: { type: 'string', description: 'path:line' },
      severity: { enum: SEVERITY },
      category: { enum: CATEGORY },
      trigger_path: { type: 'string', description: 'concrete reachable path + inputs' },
      impact: { type: 'string', description: 'observable harmful effect' },
      evidence: { type: 'string', description: 'code + reasoning' },
      confidence: { type: 'number' },
    } } } },
}

const REACHABILITY_SCHEMA = {
  type: 'object', additionalProperties: false,
  required: ['status', 'excluding_guard', 'call_path', 'explanation'],
  properties: {
    // Does any code on the calling path EXCLUDE the input/state that triggers this bug?
    status: { enum: REACH_STATUS },
    excluding_guard: { type: 'string', description: 'file:line of the check in calling code that excludes the triggering case, or "none"' },
    call_path: { type: 'string', description: 'the python -> script_interface -> core chain traced to reach this code (or where the path begins)' },
    explanation: { type: 'string', description: 'why the trigger is/ isn\'t excluded by calling code' },
  },
}

const VERDICT_SCHEMA = {
  type: 'object', additionalProperties: false,
  required: ['verdict', 'reasoning', 'reachability', 'corrected_severity'],
  properties: {
    verdict: { enum: ['confirmed', 'refuted', 'uncertain'] },
    reasoning: { type: 'string' },
    reachability: REACHABILITY_SCHEMA,
    corrected_severity: { enum: ['critical', 'high', 'medium', 'low', 'none'] },
  },
}

const SYNTH_SCHEMA = {
  type: 'object', additionalProperties: false,
  required: ['summary', 'ranked', 'systemic_patterns'],
  properties: {
    summary: { type: 'string' },
    ranked: { type: 'array', items: {
      type: 'object', additionalProperties: false,
      required: ['rank', 'title', 'file', 'severity', 'category', 'impact', 'fix_hint', 'reachability'],
      properties: {
        rank: { type: 'integer' }, title: { type: 'string' }, file: { type: 'string' },
        severity: { enum: SEVERITY }, category: { type: 'string' },
        impact: { type: 'string' }, fix_hint: { type: 'string' },
        reachability: { type: 'string', description: 'one-line: status + the excluding guard (or "none") carried from verification' },
      } } },
    systemic_patterns: { type: 'array', items: { type: 'string' } },
  },
}

// ---- prompts -------------------------------------------------------------
function scopePrompt(rootOverride) {
  return `You are SCOPING a bug sweep of the ESPResSo molecular-dynamics codebase. Do NOT review code — just enumerate.

1. Repository root: ${rootOverride
    ? 'use this absolute path: ' + rootOverride
    : 'run \`git rev-parse --show-toplevel\` and use its output'}. Call it ROOT.
2. From ROOT, enumerate EVERY C/C++ source file under src/ — extensions .cpp .hpp .cu .cuh —
   but EXCLUDE any path containing a test directory segment (unit_tests/ or tests/). Include ALL of
   src/ (src/core, src/utils, src/shapes, src/particle_observables, src/script_interface,
   src/walberla_bridge, src/scafacos, src/instrumentation, etc.) — these are ESPResSo's OWN code.
   Vendored third-party libraries live OUTSIDE src/ (e.g. libs/) and are already excluded.
   A fast way: \`cd ROOT && find src -type f \\( -name '*.cpp' -o -name '*.hpp' -o -name '*.cu' -o -name '*.cuh' \\) -not -path '*/unit_tests/*' -not -path '*/tests/*' | xargs wc -l\`
   Report each file as { path: <repo-relative, e.g. "src/core/foo.cpp">, lines: <wc -l> }.
3. Also enumerate the Python interface modules: every *.py and *.pyx under src/python (repo-relative),
   EXCLUDING test files. These drive the public-API call paths.
4. Set file_count = number of files you listed, total_lines = sum of their line counts (so the
   workflow can verify nothing was dropped). In excluded_note, say what you skipped (tests, etc.).

Return STRICTLY via the schema. Be exhaustive and exact — every source file must appear in files[].`
}

function finderPrompt(u, root) {
  const pathLines = u.paths.map(p => '- ' + root + '/' + p).join('\n')
  return `You are auditing ONE SLICE of the ESPResSo molecular-dynamics codebase for BUGS.
Repository root: ${root}

Review unit: "${u.name}"
${u.note ? 'FOCUS: ' + u.note + '\n' : ''}Files to review (this exact set — do not wander to other units):
${pathLines}

ESTABLISH GROUND TRUTH before judging anything a bug:
- Follow #includes within src/ to read the headers/contracts these files depend on.
- Read the relevant unit tests (e.g. src/core/unit_tests/, src/utils/tests/) — they encode intended behavior and prevent contract misreads.
- ESPResSo is mature and battle-tested. Assume callers respect documented preconditions UNLESS you can exhibit a real caller (Grep for callsites) that violates them.

REPORT ONLY genuine defects with a CONCRETE, REACHABLE trigger path AND an OBSERVABLE harmful impact. Bug categories:
- memory_safety (OOB, use-after-free, dangling ref/iterator, uninitialized read)
- integer_overflow / signed-unsigned / narrowing that actually misbehaves
- concurrency_mpi (wrong rank/tag, missing/incorrect reduction, race, deadlock)
- numerics (div-by-zero, NaN/Inf propagation, catastrophic precision loss ON A NORMAL PATH)
- logic (off-by-one, inverted condition, wrong operator/variable, wrong loop bound)
- resource_leak (leak, double-free, missing cleanup that accumulates)
- api_contract (invariant violated, missing error handling with real consequence)
- md_invariant (PBC/minimum-image, ghost layer, cell/Verlet list, force/energy symmetry, units)

DO NOT REPORT (these are discarded and waste reviewers' time):
- Style, naming, const-correctness, "consider adding a check", defensive nits.
- Performance — UNLESS it silently corrupts results.
- Theoretical UB/overflow with NO reachable trigger ("if N were 2^31...").
- Behavior of EXTERNAL libraries (walberla, heffte, kokkos/cabana, boost, MPI, ScaFaCoS). Review only ESPResSo's own logic and how it CALLS them.
- Anything you cannot tie to a concrete trigger + impact.

PHYSICS/NUMERICS CAUTION: Be conservative before claiming a formula is physically/mathematically wrong. Check against documented intent (comments, docs, test expected values). If unsure it is intended, either don't report it or set confidence <= 0.4.

The ask is "find ANY REAL bugs" — 5 genuine bugs beat 50 maybes. Do NOT pad. An EMPTY findings list is a perfectly good answer for a clean unit.

For each bug provide: title; file as "path:line"; severity; category; trigger_path (how execution reaches it, with what inputs); impact (what goes wrong — prefer the "silently wrong results" framing where it applies); evidence (the offending code + your reasoning); confidence 0..1 (honest).`
}

function callPathPrompt(u, root) {
  const py = u.python.map(p => '- ' + root + '/' + p).join('\n') || '  (locate via _so_name below)'
  return `You are auditing the PUBLIC-API CALL PATH of one ESPResSo subsystem: "${u.subsystem}".
Repository root: ${root}

The chain is:  Python (src/python/espressomd) -> C++ script_interface wrapper -> simulation core.
HOW THE LINK WORKS: a Python class sets \`_so_name = "Namespace::Klass"\`; the matching
script_interface class is bound via \`register_new<Klass>("Namespace::Klass")\` (grep the
subsystem's initialize.cpp). The script_interface wrapper validates/forwards parameters
(constructor, set_parameter / do_set_parameter, do_call_method) into the core object.

Layers for THIS subsystem:
Python module(s):
${py}
script_interface dir:  ${root}/${u.si}
core dir:              ${root}/${u.core}

YOUR JOB — find bugs that live IN THE SEAM, not inside a single layer:
- A parameter accepted by the Python API that is dropped, mis-typed, narrowed, or unit-mismatched
  before it reaches the core (e.g. Python passes a double the SI truncates to int; an option the
  SI silently ignores; a vector reordered/transposed in conversion).
- Validation that the Python layer DOCUMENTS or implies but NO layer actually enforces, so an
  invalid value reaches the core and corrupts the simulation or crashes.
- Conversely, a core precondition that NOTHING on the path guarantees.
- Getter/observable round-trip bugs: value set via Python that reads back wrong, or core state
  exposed to Python with a wrong sign/scale/index.
- Checkpoint/restore or MPI paths that reach the core WITHOUT going through the SI validation.

Trace at least the main parameter-setting and method-call paths end to end. Read the SI wrapper
source and the core entry points it calls. Use Grep to follow _so_name <-> register_new and to
find which core functions the SI methods invoke.

DO NOT report pure intra-core bugs (other reviewers cover those) or external-library internals.
Report ONLY genuine seam defects with a concrete trigger (a specific Python call + value) and an
observable impact. Use category "call_path" (or the most specific category that fits). An empty
list is fine if the seam is clean. For each: title; file as "path:line" (cite the layer where the
defect is); severity; category; trigger_path (the Python call + value that triggers it); impact;
evidence (the code across layers + reasoning); confidence 0..1.`
}

function verifierPrompt(f, lens, root) {
  const lensText =
    lens === 'reachability' ? 'PRIMARY LENS — REACHABILITY: Trace UP from the cited code through its callers, and where relevant all the way out to the public Python API (python -> script_interface -> core). Find a real path that reaches this code with the claimed inputs. If no reachable trigger exists, or a guard/invariant on the path excludes it, REFUTE.' :
    lens === 'correctness' ? 'PRIMARY LENS — CORRECTNESS & CONTRACT: Re-derive whether the code is actually wrong given the INTENDED contract (read headers, comments, unit tests, docs). If behavior is correct/intended, REFUTE.' :
    'Re-read the actual code at the location AND its callers and the contract (headers/tests), then decide.'
  return `You are an ADVERSARIAL verifier on the ESPResSo codebase. Repository root: ${root}
DEFAULT POSITION: this reported bug is a FALSE POSITIVE. Your job is to REFUTE it. Confirm ONLY with concrete, undeniable evidence it is real.

Reported bug:
- Title: ${f.title}
- Location: ${f.file}
- Claimed severity: ${f.severity}   Category: ${f.category}
- Trigger path: ${f.trigger_path}
- Impact: ${f.impact}
- Evidence: ${f.evidence}

${lensText}

REACHABILITY ANALYSIS (do this for EVERY finding — it is reported back to the user):
Determine whether any CALLING CODE excludes the case that triggers this bug. Trace the path the
triggering input/state would travel:
  Python defaults & arg handling (src/python/espressomd)
    -> script_interface validation (constructors, set_parameter/do_set_parameter, do_call_method;
       grep there for throw, std::invalid_argument, std::domain_error, range checks like "< 0", "must be")
    -> core-internal guards / asserts
    -> the actual callsite(s) (Grep for callers).
Decide the reachability status:
  - "unreachable": some guard on the path makes the trigger impossible via the public API and via
    other entry points (checkpoint/restore, internal callbacks, MPI). Cite the guard file:line.
  - "partially_reachable": guarded on the common API path but reachable another way (checkpoint,
    direct internal call, missing validation on one entry point).
  - "reachable": nothing on the path excludes the triggering case.
  - "unknown": you genuinely cannot determine.
Note: "guarded in one constructor" is NOT automatically "unreachable" — check checkpoint/restore,
internal callbacks, and MPI paths too.
Fill the reachability object: status, excluding_guard (the file:line of the excluding check, or
"none"), call_path (the chain you traced), explanation.

WEIGHT REACHABILITY INTO SEVERITY: if the trigger is fully unreachable, the practical severity is
much lower (latent / only reachable by future API changes or internal refactors) — reflect that in
corrected_severity. If reachable or partially reachable, keep severity driven by real-world impact.

Method: open the file at the cited location, read surrounding code, Grep for callers and the
contract, read relevant unit tests, AND inspect the calling layers above. Watch for external-library
behavior wrongly blamed on ESPResSo. Then decide:
- "confirmed": real, reachable, harmful — you can point to the exact defect and a concrete trigger.
- "refuted": not a bug (correct as written / unreachable / guaranteed by invariant / external-lib behavior / misread contract).
- "uncertain": genuinely cannot determine (treated conservatively / dropped).
Cite specific line numbers and caller evidence in reasoning. Set reachability and corrected_severity
(or "none" if refuted) AFTER the analysis.`
}

function synthPrompt(survivors) {
  return `You are the lead reviewer producing the FINAL prioritized bug list for the ESPResSo codebase (all of src/, including the python -> script_interface -> core call paths).
These candidate bugs already PASSED adversarial verification. JSON:
${JSON.stringify(survivors, null, 2)}

Produce:
1. summary: 2-4 sentences on overall code health and the most important risks.
2. ranked: EVERY bug, ranked best-first. PRIORITIZATION RULE: rank "SILENTLY corrupts simulation results / wrong physics on a normal path" ABOVE "crash/UB on rare or impossible input", even when the latter has nominally higher severity. ALSO weight REACHABILITY: a bug whose trigger is fully UNREACHABLE (excluded by calling code) is latent and should rank BELOW otherwise-comparable bugs that are reachable; reflect this in severity too. Within a tier order by confidence/agreement and blast radius. MERGE obvious duplicates (same root cause/location). Fields: rank, title, file, severity, category, impact, fix_hint (one line), and reachability (one line: carry the status + excluding guard from the input finding).
3. systemic_patterns: recurring issue classes across the codebase, if any.
Be honest — a short list is fine.`
}

// ---- partitioning (deterministic, in JS) ---------------------------------
function dirOf(p) { const i = p.lastIndexOf('/'); return i < 0 ? '' : p.slice(0, i) }

// Group manifest files by directory; split any directory whose total line count exceeds `cap`
// into stable, line-balanced chunks. Returns [{ name, kind:'files', paths:[...] }].
function partition(files, cap) {
  const byDir = new Map()
  for (const f of files) {
    const d = dirOf(f.path)
    if (!byDir.has(d)) byDir.set(d, [])
    byDir.get(d).push(f)
  }
  const units = []
  for (const dir of [...byDir.keys()].sort()) {
    const list = byDir.get(dir).slice().sort((a, b) => a.path < b.path ? -1 : 1)
    const total = list.reduce((n, f) => n + f.lines, 0)
    const label = dir.replace(/^src\//, '') || 'src'
    if (total <= cap) {
      units.push({ name: label, kind: 'files', paths: list.map(f => f.path) })
      continue
    }
    // greedy bin-pack files into chunks of <= cap lines (always at least one file per chunk)
    const chunks = []
    let cur = [], curLines = 0
    for (const f of list) {
      if (cur.length && curLines + f.lines > cap) { chunks.push(cur); cur = []; curLines = 0 }
      cur.push(f); curLines += f.lines
    }
    if (cur.length) chunks.push(cur)
    chunks.forEach((c, i) => units.push({
      name: `${label} [${i + 1}/${chunks.length}]`, kind: 'files', paths: c.map(f => f.path),
    }))
  }
  return units
}

// Build python -> script_interface -> core call-path units for every subsystem that exists in
// BOTH script_interface/ and core/. Pairs the python module of the same name when present.
function callPathUnits(files, pythonModules) {
  const dirs = new Set(files.map(f => dirOf(f.path)))
  const subOf = (prefix) => new Set([...dirs]
    .filter(d => d.startsWith(prefix) && d.slice(prefix.length).indexOf('/') < 0 && d.length > prefix.length)
    .map(d => d.slice(prefix.length)))
  const siSubs = subOf('src/script_interface/')
  const coreSubs = subOf('src/core/')
  // a few well-known python<->core name aliases (python module name : core/SI dir name)
  const alias = { electrokinetics: 'ek', analyze: 'analysis', integrate: 'integrators', thermostat: 'thermostats' }
  const pyName = (p) => p.split('/').pop().replace(/\.(py|pyx)$/, '')
  const pyBySubsystem = new Map()
  for (const p of pythonModules) {
    const n = pyName(p)
    const key = alias[n] || n
    if (!pyBySubsystem.has(key)) pyBySubsystem.set(key, [])
    pyBySubsystem.get(key).push(p)
  }
  const units = []
  for (const sub of [...siSubs].filter(s => coreSubs.has(s)).sort()) {
    units.push({
      subsystem: sub,
      si: `src/script_interface/${sub}`,
      core: `src/core/${sub}`,
      python: pyBySubsystem.get(sub) || [],
    })
  }
  return units
}

// ---- run -----------------------------------------------------------------
// Normalize args: the runtime may deliver them as a JSON string rather than an object.
let opts = args
if (typeof opts === 'string') { try { opts = JSON.parse(opts) } catch (e) { opts = {} } }
opts = opts || {}

const cap = opts.maxLinesPerUnit || 3500
const rootOverride = opts.root || null

phase('Scope')
const scope = await agent(scopePrompt(rootOverride), { label: 'scope:enumerate', phase: 'Scope', schema: SCOPE_SCHEMA })
if (!scope || !scope.files || !scope.files.length) {
  return { error: 'scope agent returned no files', scope }
}
const root = scope.root
const allFiles = scope.files
const sumLinesAll = allFiles.reduce((n, f) => n + (f.lines || 0), 0)
if (allFiles.length !== scope.file_count || sumLinesAll !== scope.total_lines) {
  log(`WARN scope self-check mismatch: listed ${allFiles.length} files / ${sumLinesAll} lines vs reported ${scope.file_count}/${scope.total_lines} — proceeding with the listed set`)
}

// Skip machine-generated source by default (e.g. walberla_bridge .../generated_kernels): it is
// auto-produced, low-signal for a human-oriented bug hunt, and a large share of the agent budget.
// Pass { includeGenerated: true } to scan it anyway.
const isGenerated = (p) => /\/generated[_-]?kernels\//.test(p) || /\bgenerated\//.test(p)
const manifest = opts.includeGenerated ? allFiles : allFiles.filter(f => !isGenerated(f.path))
const skipped = allFiles.length - manifest.length
const sumLines = manifest.reduce((n, f) => n + (f.lines || 0), 0)
if (skipped && !opts.includeGenerated) {
  const skippedLines = sumLinesAll - sumLines
  log(`Skipped ${skipped} generated-source files (${skippedLines} lines) — pass { includeGenerated: true } to scan them.`)
}

let codeUnits = partition(manifest, cap)
let cpUnits = callPathUnits(manifest, scope.python_modules || [])

// coverage assertion: every manifest file must land in exactly one code unit
const assigned = new Set(codeUnits.flatMap(u => u.paths))
const unassigned = manifest.filter(f => !assigned.has(f.path))
if (unassigned.length) log(`WARN ${unassigned.length} source files were not assigned to any unit (e.g. ${unassigned.slice(0, 3).map(f => f.path).join(', ')})`)
log(`Scoped root=${root}: ${manifest.length} src files / ${sumLines} lines -> ${codeUnits.length} code units (cap ${cap} lines) + ${cpUnits.length} call-path units. ${scope.excluded_note}`)

// optional pilot filter
if (Array.isArray(opts.only) && opts.only.length) {
  const hit = (s) => opts.only.some(tok => s.includes(tok))
  codeUnits = codeUnits.filter(u => hit(u.name))
  cpUnits = cpUnits.filter(u => hit(u.subsystem))
  log(`PILOT mode (only=${opts.only.join(',')}): ${codeUnits.length} code + ${cpUnits.length} call-path units`)
}
const maxVerifiers = opts.singleVerifier ? 1 : 2

// ---- Find + Verify (pipelined per unit) ----------------------------------
phase('Find')

// every finding (from any unit) flows through the same adversarial verify stage
const verifyFindings = (unitName, findings) => {
  if (!findings.length) return Promise.resolve({ unit: unitName, findings: [], verified: [] })
  return parallel(findings.map(f => () => {
    const lenses = (maxVerifiers === 2 && (f.severity === 'critical' || f.severity === 'high'))
      ? ['reachability', 'correctness'] : ['default']
    return parallel(lenses.map(lens => () =>
      agent(verifierPrompt(f, lens, root), { label: 'verify:' + unitName, phase: 'Verify', schema: VERDICT_SCHEMA })
    )).then(vs => ({ finding: { ...f, unit: unitName }, verdicts: vs.filter(Boolean) }))
  })).then(arr => ({ unit: unitName, findings, verified: arr.filter(Boolean) }))
}

const codeResults = pipeline(
  codeUnits,
  (u) => agent(finderPrompt(u, root), { label: 'find:' + u.name, phase: 'Find', schema: FIND_SCHEMA })
            .then(r => ({ unit: u.name, findings: (r && r.findings) || [] })),
  (res) => verifyFindings(res.unit, res.findings),
)
const cpResults = pipeline(
  cpUnits,
  (u) => agent(callPathPrompt(u, root), { label: 'callpath:' + u.subsystem, phase: 'Find', schema: FIND_SCHEMA })
            .then(r => ({ unit: 'callpath:' + u.subsystem, findings: (r && r.findings) || [] })),
  (res) => verifyFindings(res.unit, res.findings),
)
const perUnit = [...await codeResults, ...await cpResults]

// ---- decide survivors ----------------------------------------------------
const allVerified = perUnit.filter(Boolean).flatMap(r => r.verified)
const rawCount = perUnit.filter(Boolean).reduce((n, r) => n + r.findings.length, 0)
const survivors = []
for (const v of allVerified) {
  const confirmed = v.verdicts.filter(x => x.verdict === 'confirmed')
  const refuted = v.verdicts.filter(x => x.verdict === 'refuted').length
  if (confirmed.length >= 1 && confirmed.length >= refuted) {
    const sev = confirmed[0]?.corrected_severity
    // strongest exclusion claim wins (a verifier that proved unreachability outranks one that didn't)
    const order = REACH_STATUS // ['unreachable','partially_reachable','reachable','unknown']
    const reach = confirmed.map(x => x.reachability).filter(Boolean)
      .sort((a, b) => order.indexOf(a.status) - order.indexOf(b.status))[0]
      || { status: 'unknown', excluding_guard: 'none', call_path: '', explanation: '' }
    survivors.push({
      ...v.finding,
      severity: (sev && sev !== 'none') ? sev : v.finding.severity,
      reachability: reach,
      verifier_agreement: confirmed.length + '/' + v.verdicts.length,
      verifier_reasoning: v.verdicts.map(x => x.verdict + ': ' + x.reasoning),
    })
  }
}
log('raw findings: ' + rawCount + ' -> survivors after adversarial verify: ' + survivors.length)

phase('Synthesize')
let synth = { summary: 'No confirmed bugs.', ranked: [], systemic_patterns: [] }
if (survivors.length) {
  synth = await agent(synthPrompt(survivors), { label: 'synthesize', phase: 'Synthesize', schema: SYNTH_SCHEMA })
}

return {
  stats: {
    root,
    src_files: manifest.length,
    src_lines: sumLines,
    code_units: codeUnits.length,
    call_path_units: cpUnits.length,
    raw_findings: rawCount,
    survivors: survivors.length,
  },
  survivors,
  synth,
}
