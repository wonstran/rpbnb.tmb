# Memory calibration harness for rpbnb.tmb
#
# Produces the measurements that `max_workload` and the family weights in
# `rpbnb_tmb_control()` are derived from.  The published constants live in
# `TAPE_CALIBRATION` (R/utilities.R); this script is how they are obtained, and
# re-running it is how they should be checked after a template, TMB, compiler
# or allocator change.
#
# Usage
#   Rscript inst/benchmark_memory.R                     # run the grid
#   Rscript inst/benchmark_memory.R --out results.csv   # choose the output file
#   Rscript inst/benchmark_memory.R --child FAM N DRAWS THREADS   # one cell
#
# Each grid cell runs in a fresh R process so that one measurement cannot be
# contaminated by another's allocator state -- an earlier calibration was wrong
# precisely because two sweeps ran concurrently.
#
# Reported quantities, all in MiB (1 MiB = 1024^2 bytes):
#   tape_mib   retained resident-set growth across TMB::MakeADFun()
#   eval_mib   further retained growth across one obj$fn() and one obj$gr()
#   peak_mib   peak working set over the whole sequence, above baseline
# `peak_mib` governs out-of-memory failures, and it is what `max_workload` is
# calibrated against -- both the per-unit budget AND the family weights.
# `tape_mib` is reported for diagnosis only; do not derive guard constants from
# it. Frank retains 2.87x Famoye but peaks at 3.53x, and an earlier revision
# that took the weights from `tape_mib` under-budgeted Frank's peak by ~22%.

suppressMessages({
  library(rpbnb.tmb)
})

mib <- function(bytes) bytes / 1024^2

fixture <- function(family_code, n, draws) {
  set.seed(1)
  x <- rnorm(n)
  list(
    data = rpbnb.tmb:::.build_tmb_data(
      Y1 = rpois(n, 2), Y2 = rpois(n, 2),
      X1 = cbind(1, x), X2 = cbind(1, x),
      rand_idx1 = 2L, rand_idx2 = integer(0),
      Z1 = matrix(runif(draws), ncol = 1L),
      Z2 = matrix(numeric(0), draws, 0),
      dist1 = 0L, dist2 = integer(0),
      sign1 = 1L, sign2 = integer(0),
      family_code = family_code, pois1 = FALSE, pois2 = FALSE,
      lamLo = -1, lamHi = 1
    ),
    parameters = list(
      beta1 = c(0.1, 0.2), beta2 = c(-0.1, -0.15),
      log_sd1 = log(0.2), log_sd2 = numeric(0),
      log_m1 = log(0.6), log_m2 = log(0.7), z_dep = 0.1
    )
  )
}

measure_cell <- function(family_code, n, draws, threads) {
  spec <- fixture(family_code, n, draws)
  gc()
  before <- ps::ps_memory_info()
  started <- Sys.time()
  suppressWarnings(TMB::openmp(n = threads, DLL = "rpbnb.tmb"))
  obj <- TMB::MakeADFun(
    data = spec$data, parameters = spec$parameters,
    DLL = "rpbnb.tmb", silent = TRUE
  )
  after_tape <- ps::ps_memory_info()
  invisible(obj$fn(obj$par))
  invisible(obj$gr(obj$par))
  after_eval <- ps::ps_memory_info()

  data.frame(
    family_code = family_code,
    family = c("independence", "famoye", "frank", "gaussian", "clayton")[
      family_code + 2L
    ],
    n = n, draws = draws, units = as.double(n) * as.double(draws),
    threads = threads,
    tape_mib = mib(after_tape[["rss"]] - before[["rss"]]),
    eval_mib = mib(after_eval[["rss"]] - after_tape[["rss"]]),
    peak_mib = mib(after_eval[["peak_wset"]] - before[["rss"]]),
    seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
    stringsAsFactors = FALSE
  )
}

session_metadata <- function() {
  data.frame(
    r_version = paste(R.version$major, R.version$minor, sep = "."),
    platform = R.version$platform,
    tmb_version = as.character(utils::packageVersion("TMB")),
    rcppeigen_version = as.character(utils::packageVersion("RcppEigen")),
    cxx = tryCatch(
      sub("^CXX\\s*=\\s*", "", grep(
        "^CXX\\s*=", readLines(file.path(R.home("etc"), "Makeconf")),
        value = TRUE
      )[1]),
      error = function(e) NA_character_
    ),
    stringsAsFactors = FALSE
  )
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 1L && args[[1L]] == "--child") {
  cell <- measure_cell(
    as.integer(args[[2L]]), as.integer(args[[3L]]),
    as.integer(args[[4L]]), as.integer(args[[5L]])
  )
  write.csv(cell, stdout(), row.names = FALSE)
  quit(save = "no")
}

out <- if (length(args) >= 2L && args[[1L]] == "--out") {
  args[[2L]]
} else {
  "inst/extdata/memory_calibration.csv"
}

# Famoye carries the n x draws shape (cheap to tape, so the grid can be wide);
# every other family is measured at matched workloads so the family weights are
# evidence rather than extrapolation.
grid <- rbind(
  expand.grid(family_code = 0L, n = c(500L, 1000L, 2000L, 4000L),
              draws = c(50L, 100L, 200L), threads = 1L),
  expand.grid(family_code = c(-1L, 1L, 2L, 3L), n = c(1000L, 2000L),
              draws = c(100L, 200L), threads = 1L)
)

rscript <- file.path(R.home("bin"), "Rscript")
self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (!length(self)) self <- "inst/benchmark_memory.R"

results <- list()
for (i in seq_len(nrow(grid))) {
  row <- grid[i, ]
  message(sprintf(
    "[%2d/%2d] family_code=%2d n=%5d draws=%4d threads=%d",
    i, nrow(grid), row$family_code, row$n, row$draws, row$threads
  ))
  csv <- system2(
    rscript,
    c(shQuote(self), "--child", row$family_code, row$n, row$draws, row$threads),
    stdout = TRUE
  )
  results[[i]] <- utils::read.csv(text = paste(csv, collapse = "\n"))
}

results <- do.call(rbind, results)
results <- cbind(results, session_metadata())
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
utils::write.csv(results, out, row.names = FALSE)
message("wrote ", out)

# ---- Regressions the published constants come from -------------------------
# Every constant in TAPE_CALIBRATION must be printed below, marked as such, so
# a constant can be re-derived by re-running this script rather than trusted.
famoye <- subset(results, family_code == 0L)
per_units <- aggregate(cbind(tape_mib, peak_mib) ~ units, famoye, mean)
top <- utils::tail(per_units[order(per_units$units), ], 3L)

report_slope <- function(what, column) {
  model <- lm(per_units[[column]] ~ per_units$units)
  top_slope <- coef(lm(top[[column]] ~ top$units))[[2L]]
  cat(sprintf("%s_mib = %.3f + %.8g * units      R^2 = %.5f\n",
              what, coef(model)[[1L]], coef(model)[[2L]],
              summary(model)$r.squared))
  cat(sprintf("  overall slope        : %.0f bytes/unit\n",
              coef(model)[[2L]] * 1024^2))
  cat(sprintf("  large-workload slope : %.0f bytes/unit  <- published %s\n",
              top_slope * 1024^2,
              if (what == "peak") "peak_bytes_per_unit" else
                "tape_bytes_per_unit"))
  invisible(NULL)
}

cat("\n-- Famoye peak model (what max_workload budgets) ---------------\n")
report_slope("peak", "peak_mib")
cat("\n-- Famoye tape model (diagnostic only) -------------------------\n")
report_slope("tape", "tape_mib")

# Both ratios are printed because the guard was once calibrated on the wrong
# one. The PEAK column is the one family_weight comes from.
cat("\n-- Family cost relative to Famoye at matched workload ----------\n")
matched <- subset(results, units %in% subset(results, family_code != 0L)$units)
for (u in sort(unique(matched$units))) {
  base_tape <- mean(subset(matched, units == u & family_code == 0L)$tape_mib)
  base_peak <- mean(subset(matched, units == u & family_code == 0L)$peak_mib)
  for (fc in sort(unique(subset(matched, family_code != 0L)$family_code))) {
    cell <- subset(matched, units == u & family_code == fc)
    if (!nrow(cell) || !is.finite(base_peak)) next
    cat(sprintf(
      "units=%7.0f  %-13s peak %8.2f MiB  peak ratio = %.3f   (tape %.3f)\n",
      u, cell$family[1L], mean(cell$peak_mib),
      mean(cell$peak_mib) / base_peak, mean(cell$tape_mib) / base_tape
    ))
  }
}

cat("\n-- Published family_weight = ceiling(max peak ratio, 0.1) ------\n")
for (fc in sort(unique(subset(matched, family_code != 0L)$family_code))) {
  cells <- subset(matched, family_code == fc)
  ratios <- vapply(seq_len(nrow(cells)), function(i) {
    base <- mean(subset(matched, units == cells$units[i] &
                          family_code == 0L)$peak_mib)
    cells$peak_mib[i] / base
  }, numeric(1))
  cat(sprintf("%-13s max peak ratio = %.3f  ->  weight %.1f\n",
              cells$family[1L], max(ratios),
              ceiling(max(ratios) * 10) / 10))
}

cat("\n-- Peak vs retained -------------------------------------------\n")
# Name the denominator: against tape alone this ratio is far larger, and
# calling either one "the retained footprint" has misled before.
cat(sprintf("peak / (tape + eval), median: %.2f  <- published %s\n",
            median(results$peak_mib / (results$tape_mib + results$eval_mib)),
            "peak_over_tape_and_eval"))
cat(sprintf("peak / tape        , median: %.2f  (not published)\n",
            median(results$peak_mib / results$tape_mib)))
