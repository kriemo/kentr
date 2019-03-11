context('find_motifs')

fa_file <- system.file("extdata", "test.fasta", package = "kentr")

test_that("correct number of motifs are returned",{
  test_dat <- tribble(
    ~header, ~seq,
    "foo", "AATTAATTAATT",
    "bar", "GCTGNNNACNATTAANNNAJ"
  )
  query_seq <- "AATT"
  res <- find_motifs(test_dat, query_seq)

  expect_equal(nrow(res), 3)

  query_seq <- "AATT[GCT]"
  res <- find_motifs(test_dat, query_seq)
  expect_equal(nrow(res), 0)

  query_seq <- "AAT[GCTA]"
  res <- find_motifs(test_dat, query_seq)
  expect_equal(nrow(res), 3)

  query_seq <- "[A-Z]{2}ATTA"
  res <- find_motifs(test_dat, query_seq)
  expect_equal(nrow(res), 2)

})
