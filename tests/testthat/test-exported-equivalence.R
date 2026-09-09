test_that("scalar helpers stay aligned with the reference package", {
  ref <- reference_muscle_env()

  expect_equal(logg(c(-1, 0, 1, 3)), ref$logg(c(-1, 0, 1, 3)))
  expect_equal(split(11, 4), ref$split(11, 4))
  expect_equal(split(20, 7), ref$split(20, 7))

  stepf <- list(value = c(1.5, -2), left = c(1, 4), n = 6)
  expect_equal(evalStepFun(stepf), ref$evalStepFun(stepf))

  stepf_m <- list(value = matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE), left = c(1, 4), n = 6)
  expect_equal(eval_StepFun_M(stepf_m, 1), ref$eval_StepFun_M(stepf_m, 1))
  expect_equal(eval_StepFun_M(stepf_m, 2), ref$eval_StepFun_M(stepf_m, 2))

  expect_equal(Local_errors(c(3, 7), c(2, 6), 10), ref$Local_errors(c(3, 7), c(2, 6), 10))
  expect_equal(OER(c(3, 7), c(2, 6)), ref$OER(c(3, 7), c(2, 6)))
  expect_equal(FDR(c(3, 7), c(2, 6), 10), ref$FDR(c(3, 7), c(2, 6), 10))
  expect_equal(V(c(3, 7), c(2, 6), 10), ref$V(c(3, 7), c(2, 6), 10))
  expect_equal(teethfun(12, 3, h = 4), ref$teethfun(12, 3, h = 4))
})

test_that("simulation helpers match the reference package", {
  ref <- reference_muscle_env()

  expect_equal(simulQuantile_MUSCLE(32, alpha = 0.1, beta = 0.5),
               ref$simulQuantile_MUSCLE(32, alpha = 0.1, beta = 0.5))
  expect_equal(simulQuantile_MUSCLE(32, alpha = 0.3, beta = 0.5),
               ref$simulQuantile_MUSCLE(32, alpha = 0.3, beta = 0.5))
  expect_equal(simulQuantile_MUSCLE(32, alpha = 0.5, beta = 0.5),
               ref$simulQuantile_MUSCLE(32, alpha = 0.5, beta = 0.5))
  expect_equal(simulQuantile_MUSCLE(32, alpha = 0.9, beta = 0.5),
               ref$simulQuantile_MUSCLE(32, alpha = 0.9, beta = 0.5))

  set.seed(123)
  got_q <- simulQuantile_MMUSCLE(32, alpha = 0.1, beta_vec = c(0.25, 0.5, 0.75))
  set.seed(123)
  ref_q <- ref$simulQuantile_MMUSCLE(32, alpha = 0.1, beta_vec = c(0.25, 0.5, 0.75))
  expect_equal(got_q, ref_q)
})

test_that("main exported segmentation functions match the reference package", {
  ref <- reference_muscle_env()
  set.seed(42)
  y <- rnorm(64)

  q_muscle <- simulQuantile_MUSCLE(64, alpha = 0.1, beta = 0.5)
  ref_res <- ref$MUSCLE(y, q_muscle, beta = 0.5)
  got_res <- MUSCLE(y, q_muscle, beta = 0.5)
  expect_step_result_equal(ref_res, got_res)

  ref_res_split <- ref$MUSCLE(y, q_muscle, beta = 0.5, split = TRUE, m = 24)
  got_res_split <- MUSCLE(y, q_muscle, beta = 0.5, split = TRUE, m = 24)
  expect_step_result_equal(ref_res_split, got_res_split)

  ref_res_full <- ref$MUSCLE(y, q_muscle, beta = 0.5, dyadic = FALSE)
  got_res_full <- MUSCLE(y, q_muscle, beta = 0.5, dyadic = FALSE)
  expect_step_result_equal(ref_res_full, got_res_full)

  ref_res_deconv <- ref$MUSCLE(y, q_muscle, beta = 0.5, deconv = TRUE, lag = 3)
  got_res_deconv <- MUSCLE(y, q_muscle, beta = 0.5, deconv = TRUE, lag = 3)
  expect_step_result_equal(ref_res_deconv, got_res_deconv)
})

test_that("MMUSCLE matches the reference package", {
  ref <- reference_muscle_env()
  set.seed(99)
  y <- rnorm(48)
  q_matrix <- simulQuantile_MMUSCLE(48, alpha = 0.1, beta_vec = c(0.25, 0.5, 0.75))

  ref_res <- ref$MMUSCLE(y, q_matrix, beta_vec = c(0.25, 0.5, 0.75))
  got_res <- MMUSCLE(y, q_matrix, beta_vec = c(0.25, 0.5, 0.75))
  expect_step_result_equal(ref_res, got_res)
})
