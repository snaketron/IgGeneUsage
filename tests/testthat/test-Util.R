context("Tests input rules")


test_that("Null input", {


  expect_error(checkInput(usage.data = NULL,
                          mcmc.chains = NULL,
                          mcmc.cores = NULL,
                          mcmc.steps = NULL,
                          mcmc.warmup = NULL,
                          hdi.level = NULL),
               "arguments must be specified")


  expect_error(checkInput(),
               "arguments must be specified")


})



test_that("usage.data check", {
  
  expect_error(checkUsageData(usage.data = NA),
               "usage.data must be data.frame")

  expect_error(checkUsageData(usage.data = NULL),
               "usage.data must be data.frame")

  expect_error(checkUsageData(usage.data = Inf),
               "usage.data must be data.frame")

  expect_error(checkUsageData(usage.data = character()),
               "usage.data must be data.frame")

  expect_error(checkUsageData(usage.data = numeric()),
               "usage.data must be data.frame")

  expect_error(checkUsageData(usage.data = logical()),
               "usage.data must be data.frame")

  expect_error(checkUsageData(usage.data = array()),
               "usage.data must be data.frame")

  expect_error(checkUsageData(usage.data = matrix()),
               "usage.data must be data.frame")

  expect_error(checkUsageData(usage.data = data.frame()),
               "usage.data must contain the following columns: 'sample_id',
         'condition', 'gene_name' and 'gene_usage_count'")

  expect_error(checkUsageData(
    usage.data = data.frame(a = NA, b = NA, c = NA, d = NA)),
    "usage.data must contain the following columns: 'sample_id',
         'condition', 'gene_name' and 'gene_usage_count'")

  expect_error(checkUsageData(usage.data = data.frame(
    sample_id = numeric(length = 2), condition = numeric(length = 2),
    gene_name = numeric(length = 2), gene_usage_count = numeric(length = 2))),
    "column sample_id must be of character type.")

  expect_error(checkUsageData(usage.data = data.frame(
    sample_id = numeric(length = 2),
    condition = numeric(length = 2),
    gene_name = numeric(length = 2),
    gene_usage_count = numeric(length = 2),
    temp = numeric(length = 2))),
    "usage.data must contain the following columns: 'sample_id',
         'condition', 'gene_name' and 'gene_usage_count'")

  expect_error(checkUsageData(usage.data = data.frame(
    sample_id = numeric(length = 2),
    condition = numeric(length = 2),
    gene_name = numeric(length = 2))),
    "usage.data must contain the following columns: 'sample_id',
         'condition', 'gene_name' and 'gene_usage_count'")

  expect_error(checkUsageData(usage.data = data.frame(
    sample_id = numeric(length = 2),
    gene_name = numeric(length = 2))),
    "usage.data must contain the following columns: 'sample_id',
         'condition', 'gene_name' and 'gene_usage_count'")

  expect_error(checkUsageData(usage.data = data.frame(
    sample_id = numeric(length = 2))),
    "usage.data must contain the following columns: 'sample_id',
         'condition', 'gene_name' and 'gene_usage_count'")

  expect_error(checkUsageData(usage.data = data.frame(
    sample_id = numeric(length = 2),
    condition = numeric(length = 2),
    gene_name = character(length = 2),
    gene_usage_count = character(length = 2),
    stringsAsFactors = FALSE)),
    "column sample_id must be of character type.")

  expect_error(checkUsageData(usage.data = data.frame(
    sample_id = character(length = 2),
    condition = numeric(length = 2),
    gene_name = character(length = 2),
    gene_usage_count = character(length = 2),
    stringsAsFactors = FALSE)),
    "column condition must be of character type.")

  expect_error(checkUsageData(usage.data = data.frame(
    sample_id = character(length = 2),
    condition = character(length = 2),
    gene_name = character(length = 2),
    gene_usage_count = character(length = 2),
    stringsAsFactors = FALSE)),
    "column gene_usage_count must be of numeric type.")

  expect_error(checkUsageData(usage.data = data.frame(
    sample_id = character(length = 2),
    condition = character(length = 2),
    gene_name = character(length = 2),
    gene_usage_count = NA,
    stringsAsFactors = FALSE)),
    "column gene_usage_count must be of numeric type.")

  expect_error(checkUsageData(usage.data = data.frame(
    sample_id = character(length = 2),
    condition = character(length = 2),
    gene_name = character(length = 2),
    gene_usage_count = logical(length = 2),
    stringsAsFactors = FALSE)),
    "column gene_usage_count must be of numeric type.")
  
  expect_error(checkUsageData(usage.data = data.frame(
    sample_id = character(length = 2),
    condition = c("A", "A"),
    gene_name = character(length = 2),
    gene_usage_count = numeric(length = 2),
    stringsAsFactors = FALSE)),
    "exactly 2 biological conditions must be provided.")
  
  expect_error(checkUsageData(usage.data = data.frame(
    sample_id = character(length = 3),
    condition = c("A", "B", "C"),
    gene_name = character(length = 3),
    gene_usage_count = numeric(length = 3),
    stringsAsFactors = FALSE)),
    "exactly 2 biological conditions must be provided.")
})



test_that("hdi.level check", {

  expect_error(object = checkHdi(hdi.level = NA),
               regexp = "hdi\\.level must be a number in range \\(0, 1\\)")

  expect_error(object = checkHdi(hdi.level = NULL),
               regexp = "hdi\\.level must be a number in range \\(0, 1\\)")

  expect_error(object = checkHdi(hdi.level = Inf),
               regexp = "hdi\\.level must be a number in range \\(0, 1\\)")

  expect_error(object = checkHdi(hdi.level = numeric(length = 1)),
               regexp = "hdi\\.level must be a number in range \\(0, 1\\)")

  expect_error(checkHdi(hdi.level = double(length = 1)),
               regexp = "hdi\\.level must be a number in range \\(0, 1\\)")

  expect_error(checkHdi(hdi.level = integer(length = 1)),
               regexp = "hdi\\.level must be a number in range \\(0, 1\\)")

  expect_error(checkHdi(hdi.level = character(length = 1)),
               regexp = "hdi\\.level must be a number in range \\(0, 1\\)")

  expect_error(checkHdi(hdi.level = logical(length = 1)),
               regexp = "hdi\\.level must be a number in range \\(0, 1\\)")



  expect_error(checkHdi(hdi.level = numeric(length = 3)),
               regexp = "hdi\\.level must be a number in range \\(0, 1\\)")

  expect_error(checkHdi(hdi.level = double(length = 3)),
               regexp = "hdi\\.level must be a number in range \\(0, 1\\)")

  expect_error(checkHdi(hdi.level = integer(length = 3)),
               regexp = "hdi\\.level must be a number in range \\(0, 1\\)")

  expect_error(checkHdi(hdi.level = character(length = 3)),
               regexp = "hdi\\.level must be a number in range \\(0, 1\\)")

  expect_error(checkHdi(hdi.level = logical(length = 3)),
               regexp = "hdi\\.level must be a number in range \\(0, 1\\)")

})



test_that("mcmc.chains check", {

  expect_error(object = checkMcmcChains(mcmc.chains = NA),
               regexp = "mcmc.chains must be a positive integer > 0")

  expect_error(object = checkMcmcChains(mcmc.chains = NULL),
               regexp = "mcmc.chains must be a positive integer > 0")

  expect_error(object = checkMcmcChains(mcmc.chains = Inf),
               regexp = "mcmc.chains must be a positive integer > 0")

  expect_error(object = checkMcmcChains(mcmc.chains = numeric(length = 1)),
               regexp = "mcmc.chains must be a positive integer > 0")

  expect_error(object = checkMcmcChains(mcmc.chains = double(length = 1)),
               regexp = "mcmc.chains must be a positive integer > 0")

  expect_error(object = checkMcmcChains(mcmc.chains = integer(length = 1)),
               regexp = "mcmc.chains must be a positive integer > 0")

  expect_error(object = checkMcmcChains(mcmc.chains = character(length = 1)),
               regexp = "mcmc.chains must be a positive integer > 0")

  expect_error(object = checkMcmcChains(mcmc.chains = logical(length = 1)),
               regexp = "mcmc.chains must be a positive integer > 0")

  expect_silent(object = checkMcmcChains(mcmc.chains = as.integer(x = 2)))

  # len > 1
  expect_error(object = checkMcmcChains(mcmc.chains = numeric(length = 3)),
               regexp = "mcmc.chains must be a positive integer > 0")

  expect_error(object = checkMcmcChains(mcmc.chains = double(length = 3)),
               regexp = "mcmc.chains must be a positive integer > 0")

  expect_error(object = checkMcmcChains(mcmc.chains = integer(length = 3)),
               regexp = "mcmc.chains must be a positive integer > 0")

  expect_error(object = checkMcmcChains(mcmc.chains = character(length = 3)),
               regexp = "mcmc.chains must be a positive integer > 0")

  expect_error(object = checkMcmcChains(mcmc.chains = logical(length = 3)),
               regexp = "mcmc.chains must be a positive integer > 0")
})



test_that("mcmc.cores check", {

  expect_error(object = checkMcmcCores(mcmc.cores = NA),
               regexp = "mcmc.cores must be a positive integer > 0")

  expect_error(object = checkMcmcCores(mcmc.cores = NULL),
               regexp = "mcmc.cores must be a positive integer > 0")

  expect_error(object = checkMcmcCores(mcmc.cores = Inf),
               regexp = "mcmc.cores must be a positive integer > 0")

  expect_error(object = checkMcmcCores(mcmc.cores = numeric(length = 1)),
               regexp = "mcmc.cores must be a positive integer > 0")

  expect_error(object = checkMcmcCores(mcmc.cores = double(length = 1)),
               regexp = "mcmc.cores must be a positive integer > 0")

  expect_error(object = checkMcmcCores(mcmc.cores = integer(length = 1)),
               regexp = "mcmc.cores must be a positive integer > 0")

  expect_error(object = checkMcmcCores(mcmc.cores = character(length = 1)),
               regexp = "mcmc.cores must be a positive integer > 0")

  expect_error(object = checkMcmcCores(mcmc.cores = logical(length = 1)),
               regexp = "mcmc.cores must be a positive integer > 0")

  expect_silent(object = checkMcmcCores(mcmc.cores = as.integer(x = 2)))

  # len > 1
  expect_error(object = checkMcmcCores(mcmc.cores = numeric(length = 3)),
               regexp = "mcmc.cores must be a positive integer > 0")

  expect_error(object = checkMcmcCores(mcmc.cores = double(length = 3)),
               regexp = "mcmc.cores must be a positive integer > 0")

  expect_error(object = checkMcmcCores(mcmc.cores = integer(length = 3)),
               regexp = "mcmc.cores must be a positive integer > 0")

  expect_error(object = checkMcmcCores(mcmc.cores = character(length = 3)),
               regexp = "mcmc.cores must be a positive integer > 0")

  expect_error(object = checkMcmcCores(mcmc.cores = logical(length = 3)),
               regexp = "mcmc.cores must be a positive integer > 0")
})



test_that("mcmc.steps and mcmc.warmup check", {

  expect_error(object = checkMcmcSteps(mcmc.steps = NA, mcmc.warmup = NA),
               regexp = "mcmc.steps >= 500 & mcmc.warmup >= 100.")

  expect_error(object = checkMcmcSteps(mcmc.steps = NULL, mcmc.warmup = NULL),
               regexp = "mcmc.steps >= 500 & mcmc.warmup >= 100.")

  expect_error(object = checkMcmcSteps(mcmc.steps = Inf, mcmc.warmup = Inf),
               regexp = "mcmc.steps >= 500 & mcmc.warmup >= 100.")

  expect_error(object = checkMcmcSteps(mcmc.steps = numeric(length = 1),
                                       mcmc.warmup = numeric(length = 1)),
               regexp = "mcmc.steps >= 500 & mcmc.warmup >= 100.")

  expect_error(object = checkMcmcSteps(mcmc.steps = double(length = 1),
                                       mcmc.warmup = double(length = 1)),
               regexp = "mcmc.steps >= 500 & mcmc.warmup >= 100.")

  expect_error(object = checkMcmcSteps(mcmc.steps = integer(length = 1),
                                       mcmc.warmup = integer(length = 1)),
               regexp = "mcmc.steps >= 500 & mcmc.warmup >= 100.")

  expect_error(object = checkMcmcSteps(mcmc.steps = character(length = 1),
                                       mcmc.warmup = character(length = 1)),
               regexp = "mcmc.steps >= 500 & mcmc.warmup >= 100.")

  expect_error(object = checkMcmcSteps(mcmc.steps = logical(length = 1),
                                       mcmc.warmup = logical(length = 1)),
               regexp = "mcmc.steps >= 500 & mcmc.warmup >= 100.")

  expect_silent(object = checkMcmcSteps(mcmc.steps = as.integer(x = 500),
                                        mcmc.warmup = as.integer(x = 100)))

  # len > 1
  expect_error(object = checkMcmcSteps(mcmc.steps = numeric(length = 2),
                                       mcmc.warmup = numeric(length = 2)),
               regexp = "mcmc.steps >= 500 & mcmc.warmup >= 100.")

  expect_error(object = checkMcmcSteps(mcmc.steps = double(length = 2),
                                       mcmc.warmup = double(length = 2)),
               regexp = "mcmc.steps >= 500 & mcmc.warmup >= 100.")

  expect_error(object = checkMcmcSteps(mcmc.steps = integer(length = 2),
                                       mcmc.warmup = integer(length = 2)),
               regexp = "mcmc.steps >= 500 & mcmc.warmup >= 100.")

  expect_error(object = checkMcmcSteps(mcmc.steps = character(length = 2),
                                       mcmc.warmup = character(length = 2)),
               regexp = "mcmc.steps >= 500 & mcmc.warmup >= 100.")

  expect_error(object = checkMcmcSteps(mcmc.steps = logical(length = 2),
                                       mcmc.warmup = logical(length = 2)),
               regexp = "mcmc.steps >= 500 & mcmc.warmup >= 100.")
})


