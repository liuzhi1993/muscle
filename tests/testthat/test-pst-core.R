brute_window_query <- function(y, start_offset, end_offset, order,
                               lo, hi, parity) {
  starts <- seq.int(lo, hi)
  starts <- starts[starts %% 2L == parity]
  values <- vapply(starts, function(k) {
    window <- y[seq.int(k + start_offset + 1L,
                        k + end_offset + 1L)]
    sort(window, method = "radix")[[order]]
  }, numeric(1))
  max_value <- max(values)
  min_value <- min(values)
  list(
    max_value = max_value,
    max_index = starts[which(values == max_value)[1L]],
    min_value = min_value,
    min_index = starts[which(values == min_value)[1L]]
  )
}

check_pst_layout <- function(y, start_offset, end_offset, n_starts) {
  window_length <- end_offset - start_offset + 1L
  for (order in seq_len(window_length)) {
    for (lo in 0:(n_starts - 1L)) {
      for (hi in lo:(n_starts - 1L)) {
        for (parity in 0:1) {
          if (!any(seq.int(lo, hi) %% 2L == parity)) {
            next
          }
          expected <- brute_window_query(
            y, start_offset, end_offset, order, lo, hi, parity
          )
          actual <- .pst_window_query(
            y, start_offset, end_offset, n_starts,
            order, lo, hi, parity
          )
          expect_equal(actual[c("max_value", "max_index",
                                "min_value", "min_index")],
                       expected)
        }
      }
    }
  }
}

test_that("persistent range quantiles match exhaustive brute force", {
  y <- c(7, 2, 7, 4, 9, 2, 5, 9, 1, 5)
  n <- length(y)

  for (l in c(1L, 2L, 4L, 8L)) {
    check_pst_layout(y, 0L, l - 1L, n - l + 1L)

    if (l < n) {
      check_pst_layout(y, 1L, l, n - l)
    }
    if (l >= 2L && l < n) {
      check_pst_layout(y, 2L, l, n - l)
    }
  }
})

test_that("equal observations are committed in one threshold version", {
  y <- c(4, 4, 1, 4, 1, 7, 7, 1)
  result <- .pst_window_query(y, 0L, 3L, 5L, 2L, 0L, 4L, 0L)

  expect_equal(result$versions, length(unique(y)) + 1L)
  expect_equal(result[c("max_value", "max_index", "min_value", "min_index")],
               brute_window_query(y, 0L, 3L, 2L, 0L, 4L, 0L))
})
