testthat::test_that("zigzag_index()", {
  norm_dist <- sapply(seq_len(100), function(x) {
    med_i <- sample(seq_len(100), 1)
    zigzag_index(
      hist(
        rnorm(1000, med_i, med_i / sample(c(2, 3, 4, 5, 10), 1)),
        plot = FALSE
      )$counts
    )
  })
  rand_dist <- sapply(seq_len(100), function(x) {
    med_i <- sample(seq_len(100), 1)
    zigzag_index(hist(sample(seq_len(1000), 100), plot = FALSE)$counts)
  })
  testthat::expect_true(
    mean(rand_dist, na.rm = TRUE) > mean(norm_dist, na.rm = TRUE)
  )
})
