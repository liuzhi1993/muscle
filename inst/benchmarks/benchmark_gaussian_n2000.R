n <- 2000L

q_env <- new.env(parent = emptyenv())
data("q_30000_0_3_0_5", package = "musclePST", envir = q_env)
q <- q_env$q_30000_0_3_0_5[seq_len(n)]

set.seed(20260726)
mu <- c(rep(-4, 986), rep(0, 30), rep(4, 984))
y <- mu + rnorm(n, sd = sqrt(0.9))

implementations <- list(
  muscle = getExportedValue("muscle", "MUSCLE"),
  muscleAI = getExportedValue("muscleAI", "MUSCLE"),
  musclePST = getExportedValue("musclePST", "MUSCLE")
)

results <- list()
times <- vapply(names(implementations), function(package) {
  gc()
  elapsed <- system.time({
    results[[package]] <<- implementations[[package]](y, q, beta = 0.5)
  })[["elapsed"]]
  unname(elapsed)
}, numeric(1))

old <- results$muscle
comparison <- data.frame(
  package = names(times),
  elapsed = unname(times),
  speedup = unname(times[["muscle"]] / times),
  left_identical = vapply(results, function(x) identical(old$left, x$left),
                          logical(1)),
  value_equal = vapply(results, function(x) {
    isTRUE(all.equal(old$value, x$value, tolerance = 1e-12))
  }, logical(1))
)

print(comparison, row.names = FALSE)
