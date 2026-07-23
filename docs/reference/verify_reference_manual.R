manual <- "docs/reference/rpbnb.tmb-reference-manual.html"
if (!file.exists(manual)) stop("Reference-manual HTML is missing.")

html <- paste(readLines(manual, warn = FALSE), collapse = "\n")
html_code <- html
html_code <- gsub("&lt;", "<", html_code, fixed = TRUE)
html_code <- gsub("&gt;", ">", html_code, fixed = TRUE)
html_code <- gsub("&amp;", "&", html_code, fixed = TRUE)
html_code <- gsub("&quot;", "\"", html_code, fixed = TRUE)
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

ids <- unique(sub(
  '^id="|"$', "",
  regmatches(html, gregexpr('id="[^"]+"', html, perl = TRUE))[[1L]]
))
links <- unique(sub(
  '^href="#|"$', "",
  regmatches(html, gregexpr('href="#[^"]+"', html, perl = TRUE))[[1L]]
))
missing_links <- setdiff(links, ids)
if (length(missing_links)) {
  stop("Missing internal-link anchors: ", paste(missing_links, collapse = ", "))
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
  if (!grepl(snippet, html_code, fixed = TRUE)) {
    stop("Missing manual example snippet: ", snippet)
  }
}

example_blocks <- regmatches(
  html,
  gregexpr('<code class="example">[\\s\\S]*?</code>', html, perl = TRUE)
)[[1L]]
if (!length(example_blocks)) stop("No example blocks found in the manual.")

for (block in example_blocks) {
  code <- sub('^<code class="example">', "", block)
  code <- sub("</code>$", "", code)
  code <- gsub("&lt;", "<", code, fixed = TRUE)
  code <- gsub("&gt;", ">", code, fixed = TRUE)
  code <- gsub("&amp;", "&", code, fixed = TRUE)
  code <- gsub("&quot;", "\"", code, fixed = TRUE)
  tryCatch(parse(text = code), error = function(err) {
    stop("Manual example does not parse: ", conditionMessage(err))
  })
}

cat("Reference-manual structure verified.\n")
