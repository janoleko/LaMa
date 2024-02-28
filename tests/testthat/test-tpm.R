test_that("row sums are equal to one (for two states)", {
  param = rep(-1, 2)
  expect_equal(rowSums(tpm(param)), rep(1, 2))
})

test_that("row sums are equal to one (for three states)", {
  param = rep(-2, 6)
  expect_equal(rowSums(tpm(param)), rep(1, 3))
})
