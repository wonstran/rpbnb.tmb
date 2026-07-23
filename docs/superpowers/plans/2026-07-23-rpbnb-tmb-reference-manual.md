# RP-BNB TMB Reference Manual Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Create a source-controlled PDF reference manual for `rpbnb.tmb` that follows the useful structure of the CRAN maxLik reference page.

**Architecture:** Author an original, self-contained HTML manual under `docs/reference/`, using internal anchors and print CSS for a dense technical-reference layout. Generate the final PDF with headless Chrome, then render the PDF to PNG pages with Poppler for visual QA. The HTML remains the regeneration source beside the PDF.

**Tech Stack:** HTML5, CSS print media rules, Google Chrome headless PDF printing, Poppler (`pdfinfo`, `pdftoppm`), R 4.5 parsing checks, Git

## Global Constraints

- Create `docs/reference/rpbnb.tmb-reference-manual.pdf` and retain its HTML source beside it.
- Use original wording, examples, and visual styling; do not copy maxLik content, code, branding, or assets.
- Document every `NAMESPACE` export: `fit_rpbnb_tmb()`, `simulate_rpbnb_tmb()`, `rpbnb_tmb_control()`, `copula()`, `rpbnb_tmb_marginal_effects()`, and `rpbnb_tmb_elasticities()`.
- Document `print()`, `summary()`, `coef()`, `vcov()`, `logLik()`, `AIC()`, `BIC()`, and `predict()` methods for `rpbnb_tmb_fit` in one grouped entry.
- State only source-verified behavior and defaults. `n_cores` is a TMB/OpenMP thread request, not a promise that all requested threads are realized.
- State the supported dependence modes exactly: `"famoye"`, `"independence"`, and `copula("frank")`, `copula("normal")`, or `copula("kimeldorf")`.
- Use compact examples with small simulated data; do not execute expensive model fits while building the manual.
- The PDF must be readable in grayscale, have internal contents links and page numbers, and use no external images.

## File Structure

- Create: `docs/reference/rpbnb.tmb-reference-manual.html` — canonical source with content, anchors, navigation, and print CSS.
- Create: `docs/reference/rpbnb.tmb-reference-manual.pdf` — generated final artifact.
- Create: `docs/reference/verify_reference_manual.R` — source and example verification script.
- Create: `docs/reference/add_pdf_page_numbers.py` — deterministic PDF page-number stamping helper.
- Create: `tmp/pdfs/rpbnb-tmb-reference-manual/` — ignored render-QA images only; do not commit.

---

### Task 1: Create the Source Manual and Structural Verifier

**Files:**

- Create: `docs/reference/rpbnb.tmb-reference-manual.html`
- Create: `docs/reference/verify_reference_manual.R`
- Test: `docs/reference/verify_reference_manual.R`

**Interfaces:**

- Consumes: `DESCRIPTION`, `NAMESPACE`, `R/fit_rpbnb_tmb.R`, `R/simulate_rpbnb_tmb.R`, `R/marginal_effects.R`, `R/utilities.R`, and `R/methods.R`.
- Produces: an HTML document with `id` anchors for the package overview, each exported function, methods, model notes, references, and license sections; a verifier that exits nonzero on a missing export, required anchor, or non-parsing example.

- [ ] **Step 1: Write the failing structural verifier**

Create `docs/reference/verify_reference_manual.R` with this complete check:

```r
manual <- "docs/reference/rpbnb.tmb-reference-manual.html"
if (!file.exists(manual)) stop("Reference-manual HTML is missing.")

html <- paste(readLines(manual, warn = FALSE), collapse = "\n")
required_ids <- c(
  "package-overview", "fit-rpbnb-tmb", "simulate-rpbnb-tmb",
  "rpbnb-tmb-control", "copula", "marginal-effects", "elasticities",
  "fit-methods", "model-notes", "references", "license"
)
for (id in required_ids) {
  if (!grepl(paste0('id="', id, '"'), html, fixed = TRUE)) {
    stop("Missing manual anchor: ", id)
  }
}

ns <- readLines("NAMESPACE", warn = FALSE)
exports <- sub("^export\\((.*)\\)$", "\\1", grep("^export\\(", ns, value = TRUE))
for (fun in exports) {
  if (!grepl(paste0(fun, "("), html, fixed = TRUE)) {
    stop("Missing exported function in manual: ", fun)
  }
}

examples <- c(
  'sim <- simulate_rpbnb_tmb(n = 20,',
  'ctrl <- rpbnb_tmb_control(n_cores = 1L)',
  'copula("frank")'
)
for (snippet in examples) {
  if (!grepl(snippet, html, fixed = TRUE)) {
    stop("Missing manual example snippet: ", snippet)
  }
}

cat("Reference-manual structure verified.\n")
```

- [ ] **Step 2: Run the verifier and confirm the RED state**

Run:

```powershell
& 'C:\Program Files\R\R-4.5.1\bin\Rscript.exe' docs/reference/verify_reference_manual.R
```

Expected: nonzero exit with `Reference-manual HTML is missing.`

- [ ] **Step 3: Author the HTML manual**

Create `docs/reference/rpbnb.tmb-reference-manual.html` as a self-contained
HTML5 document with inline CSS. Use a document title of `rpbnb.tmb Reference
Manual`, an initial package-metadata block populated from `DESCRIPTION`, and
a contents list linking to every required `id` in the verifier.

Use the following exact function signatures in `<pre><code>` usage blocks:

```r
fit_rpbnb_tmb(formula_1, formula_2, data,
              random_1 = NULL, random_2 = NULL,
              draws = 400L, seed = 1234L, start = NULL,
              dependence = "famoye",
              control = rpbnb_tmb_control(),
              poisson_1 = FALSE, poisson_2 = FALSE)

simulate_rpbnb_tmb(n, beta1, beta2,
                   random_1 = NULL, random_2 = NULL,
                   dispersion = c(m1 = 0.5, m2 = 0.5),
                   dependence = "famoye", lambda = 0,
                   covariates = NULL, seed = NULL)

rpbnb_tmb_control(iterlim = 500L, reltol = 1e-8,
                  print_level = 0L, n_cores = 1L,
                  halton_burn = 300L)

copula(family, par = NULL)

rpbnb_tmb_marginal_effects(fit, which = c("y1", "y2", "both"),
                            type = c("AME", "MEM"), vars = NULL,
                            include_intercept = FALSE, digits = 4L, ...)

rpbnb_tmb_elasticities(fit, which = c("y1", "y2", "both"),
                        type = c("AME", "MEM"), vars = NULL,
                        include_intercept = FALSE, digits = 4L, ...)
```

For every primary reference entry, include headings in this order when they
apply: `Description`, `Usage`, `Arguments`, `Value`, `Details`, `Examples`,
and `See Also`. Include a grouped `Methods for rpbnb_tmb_fit` entry that lists
the eight registered methods and accurately notes that `predict()` currently
returns a coefficient-based prediction object rather than new-data fitted
means.

Include a compact package overview example beginning exactly with:

```r
sim <- simulate_rpbnb_tmb(n = 20,
  beta1 = c("(Intercept)" = 0.2, x1 = 0.4),
  beta2 = c("(Intercept)" = 0.1, x1 = -0.3),
  dispersion = c(m1 = 0.4, m2 = 0.5), seed = 1)
ctrl <- rpbnb_tmb_control(n_cores = 1L)
```

Mark any example requiring model fitting with `# Not run during manual build`
and do not invoke it in the verifier.

Apply these CSS requirements: US Letter `@page` margins, a sans-serif body,
monospaced code blocks with wrapping, dark text on white background, a
visible but grayscale-safe heading hierarchy, page-break avoidance inside
reference entries, a print-only footer with page counter, and `a { color:
inherit; }` for print legibility.

- [ ] **Step 4: Run the structural verifier and HTML source checks**

Run:

```powershell
& 'C:\Program Files\R\R-4.5.1\bin\Rscript.exe' docs/reference/verify_reference_manual.R
rg -n "TODO|TBD|Lorem ipsum|\[\[" docs/reference/rpbnb.tmb-reference-manual.html
```

Expected: `Reference-manual structure verified.` and no placeholder matches.

- [ ] **Step 5: Commit the source manual and verifier**

Run:

```powershell
git add -- docs/reference/rpbnb.tmb-reference-manual.html docs/reference/verify_reference_manual.R
git diff --cached --check
git commit -m "docs: add RP-BNB TMB reference manual source"
```

Expected: one commit containing only the manual source and verifier.

---

### Task 2: Generate and Visually Verify the PDF

**Files:**

- Create: `docs/reference/rpbnb.tmb-reference-manual.pdf`
- Create: `docs/reference/add_pdf_page_numbers.py`
- Verify: `docs/reference/rpbnb.tmb-reference-manual.html`
- Verify: `docs/reference/verify_reference_manual.R`

**Interfaces:**

- Consumes: the Task 1 HTML and verifier.
- Produces: a grayscale-readable, paginated PDF whose internal HTML anchors
  are represented in the contents and whose rendered pages are visually clean.

- [ ] **Step 1: Generate the PDF from the HTML source and stamp page numbers**

Create `docs/reference/add_pdf_page_numbers.py` with this complete helper:

```python
from io import BytesIO
from pathlib import Path
import sys

from pypdf import PdfReader, PdfWriter
from reportlab.pdfgen import canvas


def stamp(input_path: Path, output_path: Path) -> None:
    reader = PdfReader(str(input_path))
    writer = PdfWriter()
    total = len(reader.pages)
    for number, page in enumerate(reader.pages, start=1):
        width = float(page.mediabox.width)
        height = float(page.mediabox.height)
        packet = BytesIO()
        overlay_canvas = canvas.Canvas(packet, pagesize=(width, height))
        overlay_canvas.setFont("Helvetica", 8)
        overlay_canvas.drawCentredString(
            width / 2, 18, f"rpbnb.tmb Reference Manual - {number} of {total}"
        )
        overlay_canvas.save()
        packet.seek(0)
        page.merge_page(PdfReader(packet).pages[0])
        writer.add_page(page)
    with output_path.open("wb") as handle:
        writer.write(handle)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("usage: add_pdf_page_numbers.py INPUT.pdf OUTPUT.pdf")
    stamp(Path(sys.argv[1]), Path(sys.argv[2]))
```

Run:

```powershell
New-Item -ItemType Directory -Force tmp/pdfs/rpbnb-tmb-reference-manual | Out-Null
& 'C:\Program Files\Google\Chrome\Application\chrome.exe' --headless --disable-gpu --no-pdf-header-footer --print-to-pdf="$PWD\tmp\pdfs\rpbnb-tmb-reference-manual\raw.pdf" "$PWD\docs\reference\rpbnb.tmb-reference-manual.html"
& 'C:\Users\litabook\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' docs/reference/add_pdf_page_numbers.py tmp/pdfs/rpbnb-tmb-reference-manual/raw.pdf docs/reference/rpbnb.tmb-reference-manual.pdf
```

Expected: `docs/reference/rpbnb.tmb-reference-manual.pdf` exists and is
non-empty.

- [ ] **Step 2: Render the PDF to page images and inspect every page**

Run:

```powershell
pdfinfo docs/reference/rpbnb.tmb-reference-manual.pdf
pdftoppm -png -r 150 docs/reference/rpbnb.tmb-reference-manual.pdf tmp/pdfs/rpbnb-tmb-reference-manual/page
```

Open every `tmp/pdfs/rpbnb-tmb-reference-manual/page-*.png` at 100% zoom.
Check that the title, contents, all function entries, code blocks, section
breaks, page numbers, and final references/license pages have no clipping,
overlap, black blocks, missing glyphs, or awkward blank pages. If any issue is
found, adjust only the HTML/CSS source, regenerate the PDF, and repeat this
step.

- [ ] **Step 3: Re-run source and artifact checks**

Run:

```powershell
& 'C:\Program Files\R\R-4.5.1\bin\Rscript.exe' docs/reference/verify_reference_manual.R
pdfinfo docs/reference/rpbnb.tmb-reference-manual.pdf | Select-String -Pattern "Pages|Page size"
git diff --check
```

Expected: verifier succeeds, `pdfinfo` reports at least two pages and Letter
page dimensions, and Git reports no whitespace errors.

- [ ] **Step 4: Commit the generated PDF**

Run:

```powershell
git add -- docs/reference/add_pdf_page_numbers.py docs/reference/rpbnb.tmb-reference-manual.pdf
git diff --cached --check
git commit -m "docs: add RP-BNB TMB reference manual PDF"
```

Expected: one commit containing the PDF and its page-number helper. Do not commit
`tmp/pdfs/rpbnb-tmb-reference-manual/`.
