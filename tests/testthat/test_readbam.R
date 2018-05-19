context("read_bam")


bam_file <- system.file("extdata", "small_sorted.bam", package = "kentr")

test_that("reading bam works", {
  res <- bam_to_df(bam_file)
  expect_equal(nrow(res), 9976)
})
