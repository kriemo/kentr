context("read_tabix")

tbx_file <- system.file("extdata", "tabix.bed.gz", package = "kentr")

test_that("reading bam works", {
  res <- tabix_to_df(tbx_file)
  expect_equal(nrow(res), 9976)
})


test_that("regional query works", {
  res <- tabix_to_df(tbx_file, region = "chr1:3e6-3.1e6")
  expect_equal(nrow(res), 16)
})


test_that("missing regions return 0 row data.frame", {
  res <- tabix_to_df(tbx_file, region = "FOO")
  expect_equal(nrow(res), 0)
})

test_that("strand col is assigned",{
  res <- tabix_to_df(tbx_file)
  expect_true("strand" %in% colnames(res))
})
