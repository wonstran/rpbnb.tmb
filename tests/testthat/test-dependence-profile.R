test_that("the shared dependence link reproduces the reported natural values", {
  z <- 0.37
  lo <- -0.4
  hi <- 0.9
  coefs <- c(log_m1 = log(0.5), log_m2 = log(0.4), z_dep = z)

  for (fc in 0:3) {
    link <- .rpbnb_dependence_link(fc, lamLo = lo, lamHi = hi)
    report <- .rpbnb_natural_report(coefs, family_code = fc,
                                    lamLo = lo, lamHi = hi)
    values <- link$map(z)

    expect_identical(names(values), link$names, info = paste("family", fc))
    for (nm in link$names) {
      expect_equal(report[[nm]]$value, unname(values[[nm]]),
                   info = paste("family", fc, nm))
    }
  }
})

test_that("independence has no dependence link", {
  expect_null(.rpbnb_dependence_link(-1L))
})
