reference_muscle_root <- function() {
  testthat::test_path("..", "..", "..", "..", "muscle", "muscle")
}

reference_muscle_env <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) {
      return(cache)
    }

    root <- reference_muscle_root()
    if (!dir.exists(root)) {
      testthat::skip("reference muscle source package is not available")
    }
    env <- new.env(parent = globalenv())
    Rcpp::sourceCpp(file.path(root, "src", "MUSCLE.cpp"), env = env, verbose = FALSE)
    source(file.path(root, "R", "MUSCLE-package.R"), local = env)

    data_dir <- file.path(root, "data")
    data_files <- list.files(data_dir, pattern = "\\.rda$", full.names = TRUE)
    for (data_file in data_files) {
      load(data_file, envir = env)
    }

    cache <<- env
    cache
  }
})

expect_step_result_equal <- function(ref, got, tolerance = 1e-12) {
  expect_type(got, "list")
  expect_type(ref, "list")
  expect_equal(names(got), names(ref))
  if (!is.null(ref$n) || !is.null(got$n)) {
    expect_equal(got$n, ref$n)
  }
  if (!is.null(ref$left) || !is.null(got$left)) {
    expect_equal(got$left, ref$left)
  }
  if (!is.null(ref$value) || !is.null(got$value)) {
    expect_equal(got$value, ref$value, tolerance = tolerance)
  }
  if (!is.null(ref$first) || !is.null(got$first)) {
    expect_equal(got$first, ref$first, tolerance = tolerance)
  }
  if (!is.null(ref$last) || !is.null(got$last)) {
    expect_equal(got$last, ref$last, tolerance = tolerance)
  }
}
