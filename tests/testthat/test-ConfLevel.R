# Test that the sum of the levels add up to 1
test_that("confidence levels add up to 1", {
  expect_equal({
    alpha <- 0.1
    nd <- c(1, 2, 4, 8, 16, 32, 64, 128, 256, 512)
    ahat <- ConfLevel(nd, alpha, "bonferroni")
    sum(ahat * nd)},
    alpha)

  expect_equal({
    alpha <- 0.1
    nd <- c(1, 2, 4, 8, 16, 32, 64, 128, 256, 512)
    ahat <- ConfLevel(nd, alpha, "weighted_bonferroni")
    sum(ahat * nd)},
    alpha)
})
