# Cross-implementation checks against the reference `rpbnb` package.
#
# The data deliberately comes from THIS package's inst/extdata rather than from
# `rpbnb`'s own `rwm1984_clean`. That dataset is not schema-stable across
# releases: in rpbnb 0.2.1 it loads as a 3,874-by-1 data frame whose single
# column is named for the comma-joined header, so
#   Variable(s) not found in data: docvis, outwork, hospvis
# aborted both fits before either likelihood was compared. Because DESCRIPTION
# only Suggests `rpbnb`, that turned the package's own test result into a
# function of whether an optional dependency happened to be installed -- red on
# a machine with it, silently skipped without it.
#
# Feeding both implementations a file we control keeps this a real
# cross-implementation comparison while making the schema our own problem.
rwm1984_bnb_data <- function() {
  path <- system.file("extdata", "rwm1984_bnb.csv",
                      package = "rpbnb.tmb", mustWork = TRUE)
  d <- utils::read.csv(path)
  needed <- c("docvis", "hospvis", "outwork")
  missing <- setdiff(needed, names(d))
  if (length(missing)) {
    stop("rwm1984_bnb.csv is missing required column(s): ",
         paste(missing, collapse = ", "))
  }
  d
}

test_that("TMB Famoye logLik matches rpbnb (same data, independence)", {
  skip_on_cran()
  skip_if_not_installed("rpbnb")
  d <- rwm1984_bnb_data()

  # Independence model
  fit_rpbnb <- rpbnb::fit_bnb(
    docvis ~ outwork, hospvis ~ outwork,
    data = d, dependence = "independence"
  )
  fit_tmb <- fit_rpbnb_tmb(docvis ~ outwork, hospvis ~ outwork,
                            data = d, dependence = "independence",
                            draws = 1)  # no random coefs, draws=1 is fine

  # Compare log-likelihoods (should be very close)
  diff_ll <- abs(logLik(fit_rpbnb) - logLik(fit_tmb))
  expect_lt(diff_ll, 0.01)
})

test_that("TMB Famoye logLik matches rpbnb (famoye dependence)", {
  skip_on_cran()
  skip_if_not_installed("rpbnb")
  d <- rwm1984_bnb_data()

  fit_rpbnb <- rpbnb::fit_bnb(
    docvis ~ outwork, hospvis ~ outwork,
    data = d, dependence = "famoye",
    control = rpbnb::rpbnb_control(compute_se = FALSE)
  )
  fit_tmb <- fit_rpbnb_tmb(docvis ~ outwork, hospvis ~ outwork,
                            data = d, dependence = "famoye")

  diff_ll <- abs(logLik(fit_rpbnb) - logLik(fit_tmb))
  expect_lt(diff_ll, 0.1)  # Allow looser tolerance due to different optimizers
})
