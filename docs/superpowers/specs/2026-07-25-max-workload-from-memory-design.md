# Design: `rpbnb_tmb_max_workload()` and a memory-aware `max_workload` default

Date: 2026-07-25

## Problem

`rpbnb_tmb_control()`'s `max_workload` guard defaults to a fixed number —
`.calibration_default_workload()`, currently `700000`, derived from
`TAPE_CALIBRATION$budget_gib = 8` and `TAPE_CALIBRATION$peak_bytes_per_unit =
12083` (`R/utilities.R:31-54,99-106`). That number targets an 8 GiB peak
working set regardless of the machine actually running the fit. A user with
64 GiB available is needlessly capped at 8 GiB of budget; a user with 4 GiB
available can still pass the guard and then genuinely OOM, because the guard
was calibrated against a number that has nothing to do with their machine.

## Goal

Two additions and one behavior change:

1. `rpbnb_tmb_max_workload(budget_gib = NULL, fraction = 0.8)` — exported.
   Computes a `max_workload` value either from an explicit budget you supply,
   or from the machine's currently available memory.
2. `.detect_available_memory_gib()` — internal. Cross-platform available-memory
   detection with no new package dependency.
3. `rpbnb_tmb_control()`'s `max_workload` default changes from
   `.calibration_default_workload()` to `rpbnb_tmb_max_workload()`, so every
   fit's default budget is now machine-aware unless the caller overrides it.

Out of scope: changing `TAPE_CALIBRATION` itself, changing how the guard is
*applied* during fitting (`.check_tmb_workload()` in `R/tmb_helpers.R` is
unaffected — it still compares a computed workload against
`control$max_workload`), and auditing/pinning the ~15 existing test call
sites that rely on the default implicitly (see Risks).

## API

### `rpbnb_tmb_max_workload(budget_gib = NULL, fraction = 0.8)`

- `budget_gib` — `NULL` (default) to auto-detect, or one positive finite
  number: the memory budget in GiB you are stating explicitly.
- `fraction` — one number in `(0, 1]`, default `0.8`. Applied **only** to
  auto-detected available memory, never to an explicit `budget_gib`: a
  number you state yourself is trusted as-is, not second-guessed with a
  haircut. Ignored (with no error) when `budget_gib` is supplied.

Return value: `integer` or `double` workload units, computed as
`signif(effective_gib * 1024^3 / TAPE_CALIBRATION$peak_bytes_per_unit, 1)` —
the same formula `.calibration_default_workload()` already uses, just
parameterized by `effective_gib` instead of the fixed
`TAPE_CALIBRATION$budget_gib`.

Behavior:

| `budget_gib` | detection result | `effective_gib` | outcome |
|---|---|---|---|
| supplied | — | `budget_gib` | arithmetic on your number, no discount |
| `NULL` | succeeds (returns `g`) | `fraction * g` | arithmetic on the discounted detected value |
| `NULL` | fails (`NA_real_`) | — | **warns**, returns `.calibration_default_workload()` unchanged |

The warning on the fallback path names the reason detection is treated as
optional infrastructure, not a hard requirement: *"Could not detect
available memory on this platform; using the default 8 GiB calibration
budget. Pass `budget_gib` explicitly to set your own."*

### `.detect_available_memory_gib()`

`@keywords internal`, `@noRd`, no exported symbol. Returns one non-negative
finite double, or `NA_real_` on any failure — never errors, never blocks.

Dispatches on `Sys.info()[["sysname"]]`:

- **`"Linux"`** — reads `/proc/meminfo`, extracts the `MemAvailable:` line
  (kilobytes), converts to GiB. `MemAvailable` is the kernel's own
  reclaimable-cache-aware estimate; `MemFree` alone systematically
  undercounts memory the kernel would actually hand back on request. If
  `MemAvailable` is absent (very old kernels), fall back to `MemFree`.
- **`"Darwin"`** — runs `vm_stat`, sums `Pages free` and `Pages inactive`
  (both are safely reclaimable), multiplies by the reported page size,
  converts to GiB.
- **`"Windows"`** — runs
  `system2("wmic", c("OS", "get", "FreePhysicalMemory", "/value"),
  stdout = TRUE)`, parses the `FreePhysicalMemory=<KB>` line, converts to
  GiB.
- **Anything else, or any parse/`system2` failure on the above** —
  `NA_real_`. Every `system2()` call is wrapped in `tryCatch()`; a missing
  binary, a locale-mangled decimal separator, an empty result, or a
  sandbox that blocks subprocess execution all degrade to `NA_real_`
  rather than propagating an error into `rpbnb_tmb_control()`'s default
  argument evaluation. A default argument that can error breaks every
  caller who doesn't override it.

### `rpbnb_tmb_control()`

One-line change: `max_workload = .calibration_default_workload()` becomes
`max_workload = rpbnb_tmb_max_workload()`. Every other argument, every
validation rule, and the rest of the function body is unchanged.

## Documentation

`.calibration_doc()` (`R/utilities.R:59-97`) currently asserts the default
targets "about 8 GiB for one fit" as a flat statement, and its text is
injected into `rpbnb_tmb_control()`'s Rd via `@eval .calibration_doc()`
(`R/utilities.R:362`). That sentence becomes false once the default is
dynamic on most machines. Add one paragraph to `.calibration_doc()`'s output
explaining the new default: it auto-detects available memory and budgets 80%
of it, falling back to the fixed 8 GiB figure only when detection is
unavailable, and cross-reference `rpbnb_tmb_max_workload()`. The existing
paragraphs describing the calibration arithmetic itself (tape size, peak
ratio, family weights) are unchanged — only the "what does the default
target" claim needs to change, since that is the only part that stops being
true.

The reference manual (`docs/reference/rpbnb.tmb-reference-manual.html`),
updated in the previous session for `rpbnb_tmb_dependence_profile()`, states
`max_workload`'s default is a fixed calibrated number
(`docs/reference/rpbnb.tmb-reference-manual.html`, the `rpbnb-tmb-control`
section's `max_workload` argument entry). That needs a follow-up update once
this ships, using the same HTML/verify/regenerate-PDF pipeline established
in that session. Out of scope for this plan's tasks; recorded here so it
isn't forgotten.

## Risks, accepted explicitly

These were surfaced during design and the user approved proceeding with them
as bounded, accepted risk rather than expanding this plan's scope to close
them:

1. **Test determinism.** ~15 call sites across `tests/testthat/*.R` and
   `inst/*.R` call `rpbnb_tmb_control()` without pinning `max_workload`.
   Every one of them fits tiny models (small `n`/`draws`), so the computed
   workload will clear even an aggressively-discounted auto-detected budget
   on any realistic machine or CI runner — but the suite is no longer
   strictly deterministic across environments. Not pinned as part of this
   work.
2. **`test-parallel.R:210-225`** (`"the default is derived from the
   calibration, not restated"`) directly asserts
   `rpbnb_tmb_control()$max_workload` equals `.calibration_default_workload()`.
   That equality now only holds on the fallback (detection-failed) path.
   This test is rewritten, not merely patched to keep passing — see Testing.
3. **`test-parallel.R:227-244`** (`"the documentation is generated from the
   calibration"`) asserts `.calibration_doc()`'s text contains the formatted
   fixed default value. Once `.calibration_doc()`'s prose changes (see
   Documentation), this specific assertion is replaced; the surrounding
   assertions (peak/tape bytes-per-unit appear in the text, family weights,
   measured-family set) are unaffected and stay as-is.
4. **Detection latency.** Every call to `rpbnb_tmb_control()` without an
   explicit `max_workload` now spawns a subprocess (`wmic`/`vm_stat`) or
   reads a `/proc` file. This is expected to cost low tens of milliseconds,
   not seconds — no benchmark is planned as part of this work, but if it
   turns out to be slow in practice that is a follow-up, not a blocker here.
5. **CI/sandboxed environments.** A container with a very low memory limit
   will compute a correspondingly low `max_workload`, which could reject a
   fit that would have passed under the old fixed 8 GiB default. This is
   the *intended* behavior of a memory-aware guard (a 1 GiB container really
   shouldn't get an 8 GiB budget) but is called out because it is a genuine
   behavior change, not an oversight.

## Testing

New tests in a new file, `tests/testthat/test-max-workload.R`:

1. **Explicit `budget_gib` is pure arithmetic, no fraction applied.**
   `rpbnb_tmb_max_workload(budget_gib = 8)` equals
   `.calibration_default_workload()` exactly (same formula, same input).
   `rpbnb_tmb_max_workload(budget_gib = 16)` equals double that (within
   `signif()` rounding).
2. **Detected memory is discounted by `fraction`.** Mock
   `.detect_available_memory_gib()` via `testthat::local_mocked_bindings()`
   to return a fixed value (e.g. `10`), call
   `rpbnb_tmb_max_workload()` with the default `fraction = 0.8`, and assert
   the result equals `signif(0.8 * 10 * 1024^3 / TAPE_CALIBRATION$peak_bytes_per_unit, 1)`.
3. **`fraction` is respected when overridden.** Same mock, `fraction = 0.5`,
   assert the corresponding halved result.
4. **Failed detection warns and falls back.** Mock
   `.detect_available_memory_gib()` to return `NA_real_`, assert
   `expect_warning(rpbnb_tmb_max_workload(), "Could not detect")`, and
   assert the returned value equals `.calibration_default_workload()`.
5. **`fraction` validated.** `fraction = 0`, `fraction = 1.5`,
   `fraction = NA`, `fraction = c(0.5, 0.5)` all error.
6. **`budget_gib` validated.** Non-numeric, negative, zero, and
   non-finite `budget_gib` all error.
7. **`.detect_available_memory_gib()` never errors.** On whichever platform
   the test suite actually runs on, call it directly and assert the result
   is either a single non-negative finite double or `NA_real_` — never an
   error, never `NaN`. This is the one test that exercises the real,
   unmocked platform-detection code path, so it is inherently
   platform-dependent in *what* it returns, but not in *whether it errors*.

Rewritten in `tests/testthat/test-parallel.R`:

8. **`"the default is derived from the calibration, not restated"`
   (lines 210-225)** — replaced with two tests: one mocking
   `.detect_available_memory_gib()` to a fixed value and asserting
   `rpbnb_tmb_control()$max_workload` equals
   `rpbnb_tmb_max_workload()` computed the same way (proving the control
   default really does call through to the new function, not a stale
   copy); and one mocking detection to fail and asserting
   `rpbnb_tmb_control()$max_workload` equals `.calibration_default_workload()`
   (proving the fallback path is wired correctly end-to-end through
   `rpbnb_tmb_control()`, not just tested in isolation on
   `rpbnb_tmb_max_workload()` directly).
9. **`"the documentation is generated from the calibration"`
   (lines 227-244)** — the `expect_match(doc, format(.calibration_default_workload(), ...))`
   assertion is removed (that claim is no longer what the docs say); a new
   assertion confirms the doc text mentions `rpbnb_tmb_max_workload` and the
   word "detect", so a future edit that silently drops the dynamic-default
   explanation is caught. The other assertions in this test
   (`peak_bytes_per_unit`, `tape_bytes_per_unit`, `"frank 2.9"`,
   `measured_families` set-equality) are unchanged.

## File Structure

- Modify: `R/utilities.R` — add `.detect_available_memory_gib()` and
  `rpbnb_tmb_max_workload()` near `.calibration_default_workload()`; change
  `rpbnb_tmb_control()`'s default; extend `.calibration_doc()`'s generated
  text.
- Create: `tests/testthat/test-max-workload.R`.
- Modify: `tests/testthat/test-parallel.R` — rewrite the two tests named
  above.
- Modify: `NAMESPACE`, `man/rpbnb_tmb_max_workload.Rd` (new),
  `man/rpbnb_tmb_control.Rd` — regenerated by roxygen.
