set.seed(123)
x <- rprobMat(10, 2, 4, alpha=2)

x2_1 <- conditional(x, 2, 1)

testthat::test_that("Conditioning variable works", {
  
})