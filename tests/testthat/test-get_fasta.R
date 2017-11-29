context('get_fasta')

fa_file <- system.file("extdata", "test.fasta", package = "kentr")

test_that('basic usage works', {

  df <- data.frame(
    chrom = c("chr1",
              "chr2"),
    start = c(24000,
              49900),
    end = c(25000,
            50000))

  res <- get_sequences(df, fa_file)
  expect_true(nchar(res$seq[1]) == 1000)
  expect_true(nchar(res$seq[2]) == 100)
  expect_true(all(str_detect(res$seq[1],
                         c("a", "t", "c", "g", "A", "T", "C", "G"))))
})

test_that('strand argument works', {

  df <- data.frame(
    chrom = c("chr1",
              "chr1"),
    start = c(24000,
              24000),
    end = c(24100,
            24100),
    strand = c("+",
               "-"))

  # by default with strand col, will return rev comp for -
  res <- get_sequences(df, fa_file)
  expect_true(res$seq[1] == revComp(res$seq[2]))

  res <- get_sequences(df, fa_file, strand = F)
  expect_true(res$seq[1] == res$seq[2])

})

