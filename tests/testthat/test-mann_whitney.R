context("mann_whitney")

test_normal <- sapply(1:100, function(x) {
  set.seed(x)
  ids <- rep(c("A", "B"), each = 50 + x)
  nums1 <- rnorm(50 + x)
  set.seed(x + 1)
  nums2 <- nums1 + runif(1, -1, 1)
  df <- data.frame(ids = ids,
                   num = c(nums1, nums2))

  df$num2 <- df$num
  df$num3 <- df$num

  a <- mann_whitney(df, "ids")
  b <- lapply(c("num",
                "num2",
                "num3"),
              function(x) {
                a <- wilcox.test(df[df$ids == "A", x],
                                 df[df$ids == "B", x])
                a$p.value})
  a$pval == b
})

test_that("stats::wilcox.test == mann_whitney for n > 50", {
  # uses normal approximation at > 50 observations
  expect_true(all(test_normal))
})

test_exact <- sapply(1:50, function(x) {
  set.seed(x)
  ids <- rep(c("A", "B"), each = x)
  nums1 <- rnorm(x)
  set.seed(x + 1)
  nums2 <- nums1 + runif(1, -1, 1)
  df <- data.frame(ids = ids,
                   num = c(nums1, nums2))

  df$num2 <- df$num
  df$num3 <- df$num

  a <- mann_whitney(df, "ids")
  b <- lapply(c("num",
                "num2",
                "num3"),
              function(x) {
                a <- wilcox.test(df[df$ids == "A", x],
                                 df[df$ids == "B", x])
                a$p.value})
  a$pval == b
})

test_that("stats::wilcox.test == mann_whitney for n <= 50", {
  # produces exact p-values
  expect_true(all(test_exact))
})

test_ties <- sapply(seq(4, 100, 2),
                    function(x) {
  set.seed(x)
  ids <- rep(c("A", "B"), each = x)
  nums1 <- rnorm(as.integer(x / 2))
  nums1 <- c(nums1, nums1)
  set.seed(x + 1)
  nums2 <- nums1 + runif(1, -1, 1)
  df <- data.frame(ids = ids,
                   num = c(nums1, nums2))

  df$num2 <- df$num
  df$num3 <- df$num

  ## will warn when n < 50 and ties are present
  a <- suppressWarnings(mann_whitney(df, "ids"))
  b <- lapply(c("num",
                "num2",
                "num3"),
              function(x) {
                a <- suppressWarnings(wilcox.test(df[df$ids == "A", x],
                                 df[df$ids == "B", x]))
                a$p.value})
  a$pval == b
})

test_that("stats::wilcox.test == mann_whitney in the presence of ties", {
  # produces exact p-values
  expect_true(all(test_ties))
})

test_that("warn when n < 50 and ties", {
  ids <- rep(c("A", "B"), each = 10)
  nums <- rep(1:10, each = 2)
  df <- data.frame(ids, nums)
  expect_warning(mann_whitney(df, "ids"))
})

test_that("NAs are handled", {
  ids <- rep(c("A", "B"), each = 10)
  nums1 <- c(1:5, rep(NA, 5))
  nums2 <- nums1 + 1
  df <- data.frame(ids, nums = c(nums1, nums2))
  expect_warning(mw_test <- mann_whitney(df, "ids"))
  suppressWarnings(wtest <- wilcox.test(nums ~ ids, df)$p.value)
  expect_equal(mw_test$pval, wtest)
})

test_that("infinite values are dropped", {
  ids <- rep(c("A", "B"), each = 10)
  nums1 <- c(1:5, rep(Inf, 5))
  nums2 <- nums1 + 1
  df <- data.frame(ids, nums = c(nums1, nums2))
  expect_warning(mw_test <- mann_whitney(df, "ids"))
  suppressWarnings(wtest <- wilcox.test(nums ~ ids, df)$p.value)
  expect_equal(mw_test$pval, wtest)
})

test_that("all Inf vectors are caught", {
  ids <- rep(c("A", "B"), each = 10)
  nums1 <- rep(Inf, 10)
  nums2 <- nums1 + 1
  df <- data.frame(ids,
                   nums = c(nums1, nums2),
                   nums2 = rnorm(10))
  expect_warning(mw_test <- mann_whitney(df, "ids"))
  expect_true(is.na(mw_test$pval[1]))
  expect_true(is.na(mw_test$w_stat[1]))
})
