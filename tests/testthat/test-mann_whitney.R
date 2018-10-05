context("mann_whitney")

test_res <- sapply(1:100, function(x) {
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
  expect_true(all(test_res))
})
