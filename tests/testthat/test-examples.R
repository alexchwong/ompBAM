source("wrap_examples.R")

.test_ParaBAM <- function() {
  expect_equal(test_example(), 0)
}
test_that("test_ParaBAM", {
  .test_ParaBAM()
})