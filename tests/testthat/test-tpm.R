test_that("row sums are equal to one (for two states)", {
  param = rep(-1, 2)
  expect_equal(rowSums(tpm(param)), c(S1 = 1, S2 = 1))
})

test_that("row sums are equal to one (for three states)", {
  param = rep(-2, 6)
  expect_equal(rowSums(tpm(param)), c(S1 = 1, S2 = 1, S3 = 1))
})
