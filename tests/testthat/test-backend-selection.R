test_that("MUSCLE selects each backend and defaults to PRT", {
  y <- c(rep(0, 4), rep(3, 4), rep(0, 4), rep(-2, 4))
  q <- rep(0, length(y))

  default <- MUSCLE(y, q)
  prt <- MUSCLE(y, q, implementation = "PRT")
  pst <- MUSCLE(y, q, implementation = "PST")
  original <- MUSCLE(y, q, implementation = "original")

  expect_equal(default, prt)
  expect_equal(pst$left, original$left)
  expect_equal(pst$value, original$value)
})

test_that("invalid backend names are rejected", {
  expect_error(
    MUSCLE(rep(0, 4), rep(0, 4), implementation = "unknown"),
    "implementation"
  )
  expect_error(
    MUSCLE(rep(0, 4), rep(0, 4), implementation = "AI"),
    "implementation"
  )
  expect_error(
    MUSCLE(rep(0, 4), rep(0, 4), implementation = "MUSCLE"),
    "implementation"
  )
})

test_that("cached MMUSCLE backends match original", {
  set.seed(42)
  y <- rnorm(48)
  beta_vec <- c(0.25, 0.5, 0.75)
  q_matrix <- simulQuantile_MMUSCLE(48, alpha = 0.1, beta_vec = beta_vec)

  reference <- MMUSCLE(y, q_matrix, beta_vec = beta_vec,
                       implementation = "original")
  for (implementation in c("PRT", "PST")) {
    result <- MMUSCLE(y, q_matrix, beta_vec = beta_vec,
                      implementation = implementation)
    expect_equal(result$left, reference$left)
    expect_equal(result$value, reference$value)
    expect_equal(result$first, reference$first)
    expect_equal(result$last, reference$last)
  }
})

test_that("single-beta MMUSCLE agrees with MUSCLE", {
  set.seed(123)
  y <- rnorm(48)
  beta <- 0.5
  q <- simulQuantile_MUSCLE(length(y), alpha = 0.1, beta = beta)
  q_matrix <- matrix(q, nrow = 1L)

  for (implementation in c("PRT", "PST", "original")) {
    muscle_result <- MUSCLE(y, q, beta = beta,
                            implementation = implementation)
    mmuscle_result <- MMUSCLE(y, q_matrix, beta_vec = beta,
                              implementation = implementation)
    expect_equal(mmuscle_result$left, muscle_result$left)
    expect_equal(as.numeric(mmuscle_result$value[1, ]), muscle_result$value)
    expect_equal(as.numeric(mmuscle_result$first[1, ]), muscle_result$first)
    expect_equal(as.numeric(mmuscle_result$last[1, ]), muscle_result$last)
  }
})
