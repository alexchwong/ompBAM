.test_idxstats <- function(threads, dataset) {
    require(ompBAMExample)
    idxstats <- getFromNamespace("idxstats_pbam", "ompBAMExample")
    return(idxstats(example_BAM(dataset),threads, TRUE))
}

.test_ompBAM <- function() {
  expect_equal(.test_idxstats(1, "Unsorted"), 0)
  expect_equal(.test_idxstats(2, "scRNAseq"), 0)
}

test_that("test_ompBAM", {
  install_ompBAM_example()
  .test_ompBAM()
})
