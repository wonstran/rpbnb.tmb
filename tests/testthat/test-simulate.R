simulation_coefficients <- function() {
  list(
    beta1 = c("(Intercept)" = 0.2, x = 0.3),
    beta2 = c("(Intercept)" = 0.1, x = -0.2)
  )
}

test_that("copula simulation ignores the Famoye-only lambda argument", {
  beta <- simulation_coefficients()
  dep <- copula("kimeldorf", par = 5)

  default_lambda <- simulate_rpbnb_tmb(
    n = 300,
    beta1 = beta$beta1,
    beta2 = beta$beta2,
    dependence = dep,
    seed = 2718
  )
  explicit_lambda <- simulate_rpbnb_tmb(
    n = 300,
    beta1 = beta$beta1,
    beta2 = beta$beta2,
    dependence = dep,
    lambda = 1,
    seed = 2718
  )

  expect_identical(default_lambda$data, explicit_lambda$data)
  expect_gt(cor(default_lambda$data$y1, default_lambda$data$y2,
                method = "kendall"), 0.25)
})

test_that("simulation rejects unrecognized dependence values directly", {
  beta <- simulation_coefficients()

  expect_error(
    simulate_rpbnb_tmb(
      n = 10,
      beta1 = beta$beta1,
      beta2 = beta$beta2,
      dependence = "garbage",
      seed = 1
    ),
    '`dependence` must be "famoye", "independence", or copula\\(\\)'
  )
})

test_that("Gaussian copula simulation reports its optional dependency", {
  beta <- simulation_coefficients()

  if (!requireNamespace("pbivnorm", quietly = TRUE)) {
    expect_error(
      simulate_rpbnb_tmb(
        n = 10,
        beta1 = beta$beta1,
        beta2 = beta$beta2,
        dependence = copula("normal", par = 0.5),
        seed = 1
      ),
      "pbivnorm"
    )
  } else {
    expect_no_error(
      simulate_rpbnb_tmb(
        n = 10,
        beta1 = beta$beta1,
        beta2 = beta$beta2,
        dependence = copula("normal", par = 0.5),
        seed = 1
      )
    )
  }
})
