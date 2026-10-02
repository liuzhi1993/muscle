test_that("dyadic MUSCLE matches the reference across signals and quantiles", {
  ref <- reference_muscle_env()
  signals <- list(
    constant = rep(3, 24),
    ties = rep(c(-2, -2, 0, 0, 4, 4), 4),
    piecewise = c(rep(-4, 10), rep(0, 8), rep(4, 12)),
    outliers = c(20, rep(0, 14), -20, rep(1, 14), 30, -30)
  )

  for (beta in c(0.25, 0.5)) {
    for (y in signals) {
      q <- simulQuantile_MUSCLE(length(y), alpha = 0.3, beta = beta)
      expect_step_result_equal(
        ref$MUSCLE(y, q, beta = beta),
        MUSCLE(y, q, beta = beta)
      )
    }
  }
})

test_that("high-quantile shifted boundary is deterministic", {
  signals <- list(
    rep(c(-2, 0, 4), 8),
    c(rep(-4, 10), rep(0, 8), rep(4, 12))
  )

  for (beta in c(0.75, 0.9)) {
    for (y in signals) {
      q <- simulQuantile_MUSCLE(length(y), alpha = 0.3, beta = beta)
      result <- MUSCLE(y, q, beta = beta)
      expect_type(result, "list")
      expect_named(result, c("value", "left", "n", "first", "last"))
      expect_equal(result$n, length(y))
    }
  }
})

test_that("dyadic MUSCLE matches the reference on randomized inputs", {
  ref <- reference_muscle_env()

  for (seed in 1:12) {
    set.seed(seed)
    n <- c(16L, 24L, 32L)[(seed - 1L) %% 3L + 1L]
    y <- rnorm(n)
    if (seed %% 3L == 0L) {
      y <- round(y, 1)
    }
    beta <- c(0.25, 0.5)[(seed - 1L) %% 2L + 1L]
    q <- simulQuantile_MUSCLE(n, alpha = 0.3, beta = beta)

    expect_step_result_equal(
      ref$MUSCLE(y, q, beta = beta),
      MUSCLE(y, q, beta = beta)
    )
  }
})

test_that("dyadic deconvolution path matches the reference", {
  ref <- reference_muscle_env()

  for (seed in 21:24) {
    set.seed(seed)
    y <- round(rnorm(32), 2)
    beta <- c(0.25, 0.5, 0.75, 0.5)[seed - 20L]
    q <- simulQuantile_MUSCLE(32, alpha = 0.3, beta = beta)

    expect_step_result_equal(
      ref$MUSCLE(y, q, beta = beta, deconv = TRUE, lag = 3),
      MUSCLE(y, q, beta = beta, deconv = TRUE, lag = 3)
    )
  }
})
