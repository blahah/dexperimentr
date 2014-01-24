context('Conditions')

test_that("all_conditions generates all possible sets of conditions", {
  samples <- c(1,1,2,2,3,3)
  possible_conditions <- list(
    c(list(1, 2, 3)),
    c(list(1,2), list(3)),
    c(list(1), list(2,3)),
    c(list(1), list(2), list(3))
  )
  expect_that(get_package('datasets'), is_true(), label="failed to load datasets package")
  expect_that(get_package('stats'), is_true(), label="failed to load stats package")
})