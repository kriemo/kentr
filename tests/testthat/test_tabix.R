context("read_tabix")

tbx_file <- system.file("extdata", "tabix.bed.gz", package = "kentr")

test_that("reading tabix works", {
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

test_that("chroms can be listed", {
  res <- get_tabix_chroms(tbx_file)
  expect_true(all(res == c("chr1", "chr2")))
})
