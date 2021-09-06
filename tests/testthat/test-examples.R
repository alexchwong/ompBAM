source("wrap_examples.R")

.test_ompBAM <- function() {
  expect_equal(test_example(), 0)
}
test_that("test_ompBAM", {
  .test_ompBAM()
})