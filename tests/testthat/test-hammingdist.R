context('hamming_dist')


test_that('hamming distance count is correct', {

  a = "AAAAT"
  b = "AAAAA"
  res <- hamming_distance_iterator(a, c(a, b))
  expect_true(all(res == c(0, 1)))

})

test_that('mismatch positions are reported correctly', {

  a = "AAAAT"
  b = "AAAAA"
  res <- mismatch_position_calc(a, c(a, b))
  expect_true(all(res == c(0, 0, 0, 0, 1)))

})
