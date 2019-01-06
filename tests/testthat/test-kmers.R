context('kmers')

fa_file <- system.file("extdata", "test.fasta", package = "kentr")

df <- data.frame(
  chrom = c("chr1",
            "chr1"),
  start = c(24000,
            24000),
  end = c(24100,
          24100),
  strand = c("+",
             "-"))

res <- get_sequences(df, fa_file)

test_that('basic usage works', {

  kmers <- get_kmers(res$seq, 3)
  expect_true(is.list(kmers))
  # all length 3
  expect_true(all(sapply(kmers[[1]]$kmer, nchar) == 3))
  # all uppercase
  expect_true(all( toupper(kmers[[1]]$kmer) == kmers[[1]]$kmer))
})


test_that('homopolymers are correct', {
  test_seq <- c("AAAAAA")
  kmers <- get_kmers(test_seq, 1)

  expect_true(kmers[[1]]$kmer == "A")
  expect_true(kmers[[1]]$counts == 6)

  test_seq <- c("AAAAAA")
  kmers <- get_kmers(test_seq, 2)
  expect_true(kmers[[1]]$kmer == "AA")
  expect_true(kmers[[1]]$counts == 5)

  test_seq <- c("AAAAAA")
  kmers <- get_kmers(test_seq, 6)
  expect_true(kmers[[1]]$kmer == "AAAAAA")
  expect_true(kmers[[1]]$counts == 1)

})

test_that('NA reported for seq too small', {
  test_seq <- c("A")
  kmers <- get_kmers(test_seq, 10)

  expect_true(all(is.na(kmers[[1]])))

})



test_that('both_strands arg works', {
  test_seq <- c("AATTAA")
  kmers <- get_kmers(test_seq, 6)

  rc_test_seq <- revComp(test_seq)
  rc_kmers <- get_kmers(rc_test_seq, 6, both_strands = TRUE)

  expect_equal(kmers, rc_kmers)

  test_seq <- c("GGAACCTT")
  kmers <- get_kmers(test_seq, 2)
  expect_true("GG" %in% kmers[[1]]$kmer)
  expect_true("TT" %in% kmers[[1]]$kmer)

  kmers <- get_kmers(test_seq, 2, both_strands = TRUE)
  expect_false("GG" %in% kmers[[1]]$kmer)
  expect_false("TT" %in% kmers[[1]]$kmer)

  kmer_df <- kmers[[1]]
  expect_true(kmer_df[kmer_df$kmer == "AA", "counts"] == 2)
  expect_true(kmer_df[kmer_df$kmer == "CC", "counts"] == 2)
})
