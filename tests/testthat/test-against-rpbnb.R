test_that("TMB Famoye logLik matches rpbnb (same data, independence)", {
  skip_on_cran()
  skip_if_not_installed("rpbnb")
  library(rpbnb)

  data(rwm1984_clean, package = "rpbnb")
  d <- rwm1984_clean

  # Independence model
  fit_rpbnb <- fit_bnb(docvis ~ outwork, hospvis ~ outwork,
                       data = d, dependence = "independence")
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
  library(rpbnb)

  data(rwm1984_clean, package = "rpbnb")
  d <- rwm1984_clean

  fit_rpbnb <- fit_bnb(docvis ~ outwork, hospvis ~ outwork,
                       data = d, dependence = "famoye",
                       control = rpbnb_control(compute_se = FALSE))
  fit_tmb <- fit_rpbnb_tmb(docvis ~ outwork, hospvis ~ outwork,
                            data = d, dependence = "famoye")

  diff_ll <- abs(logLik(fit_rpbnb) - logLik(fit_tmb))
  expect_lt(diff_ll, 0.1)  # Allow looser tolerance due to different optimizers
})
