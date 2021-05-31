source("wrap_examples.R")

.test_ParaBAM <- function() {
  expect_error(test_example(), NA)
}
test_that("test_ParaBAM", {
  .test_ParaBAM()
})