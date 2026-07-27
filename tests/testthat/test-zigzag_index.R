testthat::test_that("zigzag_index()", {
  norm_dist <- c(1, 2, 3, 5, 3, 2, 1) |> zigzag_index()
  rand_dist <- c(1, 6, 5, 3, 9, 5, 4) |> zigzag_index()
  message(c("rand_dist: ", rand_dist))
  message(c("norm_dist: ", norm_dist))
  testthat::expect_true(rand_dist > norm_dist)
})

testthat::test_that("zigzag_index() < 4", {
  norm_dist <- c(1, 1, 2) |> zigzag_index()
  testthat::expect_true(is.na(norm_dist))
})
