testthat::test_that("zigzag_index()", {
  norm_dist <- c(1, 1, 2, 3, 5, 10, 5, 8, 2, 1, 1) |> zigzag_index()
  rand_dist <- c(1, 6, 5, 3, 9, 5, 6, 1, 8, 3, 5, 8) |> zigzag_index()
  testthat::expect_true(
    (rand_dist[is.finite(rand_dist)] |> mean()) >
      (norm_dist[is.finite(norm_dist)] |> mean())
  )
})
