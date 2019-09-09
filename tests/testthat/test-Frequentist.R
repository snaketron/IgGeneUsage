context("Tests Frequentist methods")


test_that("T-test check", {

  dummy.out <- data.frame(t.test.pvalue = NA,
                          t.test.tvalue = NA,
                          t.test.L95 = NA,
                          t.test.H95 = NA,
                          t.test.fdr.pvalue = as.numeric(NA),
                          stringsAsFactors = FALSE)

  # CS: 0
  usage.data <- list(Y = matrix(data = numeric(length = 1), nrow = 1, ncol = 1),
                     X = c("a"), N = c(1), N_gene = 1)
  out <- getTTestStats(usage.data = usage.data)
  expect_equal(object = out, expected = dummy.out)


  # CS: 1
  usage.data <- list(Y = matrix(data = numeric(length = 2), nrow = 1, ncol = 2),
                     X = c("a", "b"), N = c(1, 1), N_gene = 1)
  out <- getTTestStats(usage.data = usage.data)
  expect_equal(object = out, expected = dummy.out)


  # CS: 2
  usage.data <- list(Y = matrix(data = numeric(length = 3), nrow = 1, ncol = 3),
                     X = c("a", "b", "b"), N = c(1, 1, 1), N_gene = 1)
  out <- getTTestStats(usage.data = usage.data)
  expect_equal(object = out, expected = dummy.out)


  dummy.out <- data.frame(t.test.pvalue = NaN,
                          t.test.tvalue = NaN,
                          t.test.L95 = NaN,
                          t.test.H95 = NaN,
                          t.test.fdr.pvalue = NaN,
                          stringsAsFactors = FALSE)
  rownames(dummy.out) <- "t"


  # CS: 3
  usage.data <- list(Y = matrix(data = numeric(length = 4), nrow = 1, ncol = 4),
                     X = c("a", "a", "b", "b"), N = c(1, 1, 1, 1), N_gene = 1)
  out <- getTTestStats(usage.data = usage.data)
  expect_equal(object = out, expected = dummy.out)
})






test_that("U-test check", {


  # CS: 0
  dummy.out <- data.frame(u.test.pvalue = NA,
                          u.test.wvalue = NA,
                          u.test.fdr.pvalue = as.numeric(NA),
                          stringsAsFactors = FALSE)
  usage.data <- list(Y = matrix(data = numeric(length = 1), nrow = 1, ncol = 1),
                     X = c("a"), N = c(1), N_gene = 1)
  out <- getManUStats(usage.data = usage.data)
  expect_equal(object = out, expected = dummy.out)


  # CS: 1
  dummy.out <- data.frame(u.test.pvalue = NaN,
                          u.test.wvalue = 0.5,
                          u.test.fdr.pvalue = NaN,
                          stringsAsFactors = FALSE)
  rownames(dummy.out) <- "W"
  usage.data <- list(Y = matrix(data = numeric(length = 2), nrow = 1, ncol = 2),
                     X = c("a", "b"), N = c(1, 1), N_gene = 1)
  out <- getManUStats(usage.data = usage.data)
  expect_equal(object = out, expected = dummy.out)


  # CS: 2
  dummy.out <- data.frame(u.test.pvalue = NaN,
                          u.test.wvalue = 1,
                          u.test.fdr.pvalue = NaN,
                          stringsAsFactors = FALSE)
  rownames(dummy.out) <- "W"
  usage.data <- list(Y = matrix(data = numeric(length = 3), nrow = 1, ncol = 3),
                     X = c("a", "b", "b"), N = c(1, 1, 1), N_gene = 1)
  out <- getManUStats(usage.data = usage.data)
  expect_equal(object = out, expected = dummy.out)
})

