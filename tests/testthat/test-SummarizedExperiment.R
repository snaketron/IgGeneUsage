context("Tests Summarized Experiment check")


test_that("Summarized Experiment check", {
  # CS: 0
  counts <- matrix(data = 0, nrow = 6, ncol = 6)
  colnames(counts) <- LETTERS[1:6]
  rownames(counts) <- as.character(1:6)
  se.0 <- SummarizedExperiment(assays = list(counts = counts))
  expect_error(convertSummarizedExperiment(usage.data.se = se.0),
               regexp = "colData\\(usage.data\\) is empty")
  
  
  # CS: 1
  counts <- matrix(data = 0, nrow = 6, ncol = 6)
  colnames(counts) <- LETTERS[1:6]
  rownames(counts) <- as.character(1:6)
  colData <- DataFrame(condition = rep(c("A", "B"), each = 3),
                       sample_id = LETTERS[1:6],
                       row.names = LETTERS[1:6])
  se.1 <- SummarizedExperiment(assays = list(counts = counts), 
                               colData = colData)
  expect_silent(object = convertSummarizedExperiment(usage.data.se = se.1))
  expect_is(object = convertSummarizedExperiment(usage.data.se = se.1), 
            class = "data.frame")
})
