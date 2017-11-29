context("read_bam")


bam_file <- system.file("extdata", "small_sorted.bam", package = "kentr")

test_that("reading bam works", {
  res <- read_bam(bam_file)
  expect_equal(nrow(res), 9976)
})
