context('Dependencies')

test_that("get_package loads installed packages", {
  expect_that(get_package('datasets'), is_true(), label="failed to load datasets package")
  expect_that(get_package('stats'), is_true(), label="failed to load stats package")
})

test_that("get_package installs missing packages"), {
  expect_that(get_package('EBSeq', bioconductor=TRUE))
}