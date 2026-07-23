# RP-BNB TMB Reference Manual Design

## Goal

Create a standalone PDF reference manual for `rpbnb.tmb`, organized after the
CRAN `maxLik` HTML reference manual. The manual is a package reference, not a
tutorial or a verbatim copy of the `maxLik` page.

## Deliverable

Create `docs/reference/rpbnb.tmb-reference-manual.pdf`. Its structure follows
the useful parts of the `maxLik` template:

1. package title, version, author, maintainer, description, and dependencies;
2. a linked contents page;
3. a package overview with a short, runnable fitting example;
4. grouped function reference entries;
5. a methods reference entry;
6. model, dependence, and random-coefficient notes;
7. package citation and licensing information.

The PDF must use original wording and examples derived from this repository.
It must not copy text, code, branding, or visual assets from `maxLik`.

## Reference Entries

The manual documents all exported functions from `NAMESPACE`:

- `fit_rpbnb_tmb()`;
- `simulate_rpbnb_tmb()`;
- `rpbnb_tmb_control()`;
- `copula()`;
- `rpbnb_tmb_marginal_effects()`;
- `rpbnb_tmb_elasticities()`.

It also documents the methods registered for `rpbnb_tmb_fit`: `print()`,
`summary()`, `coef()`, `vcov()`, `logLik()`, `AIC()`, `BIC()`, and `predict()`.

Each primary function entry uses the same predictable sequence as the CRAN
reference page where applicable: Description, Usage, Arguments, Value,
Details, Examples, and See Also. Method documentation is grouped in one entry
to avoid repetition.

## Content Sources and Boundaries

Function signatures, defaults, object fields, and behavior come from the R
source files and `DESCRIPTION` in this repository. Examples are compact and
use small simulated data; expensive estimation examples are marked as examples
only and are not executed during document generation.

The manual states the supported dependence modes accurately: `"famoye"`,
`"independence"`, and `copula("frank")`, `copula("normal")`, or
`copula("kimeldorf")`. It identifies `n_cores` as the TMB/OpenMP thread
request and does not promise that every requested core will be realized.

## PDF Design and Navigation

Use a dense technical-reference style: restrained typography, monospaced usage
and example blocks, a contents page, page numbers, consistent entry headings,
and internal PDF links from the contents to entries. The reference is readable
when printed in grayscale and uses no external images.

## Verification

Before delivery:

- check every exported function and registered method against `NAMESPACE`;
- verify all displayed function signatures and defaults against source;
- verify every example parses as R code without running expensive fits;
- render the PDF and inspect every page for clipping, missing characters,
  broken code blocks, page-number problems, and contents links;
- leave source artifacts alongside the PDF so the manual can be regenerated.
