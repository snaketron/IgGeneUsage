context("Tests input rules")


test_that("Null input", {


  expect_error(checkInput(usage.data = NULL,
                          mcmc.chains = NULL,
                          mcmc.cores = NULL,
                          mcmc.steps = NULL,
                          mcmc.warmup = NULL,
                          hdi.level = NULL),
               "arguments must be specified")


  expect_error(checkInput(usage.data = NA,
                          mcmc.chains = NA,
                          mcmc.cores = NA,
                          mcmc.steps = NA,
                          mcmc.warmup = NA,
                          hdi.level = NA),
               "arguments must be specified")


  expect_error(checkInput(),
               "arguments must be specified")


})



test_that("usage.data", {

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

  expect_error(checkUsageData(usage.data = data.frame(
    a = NA, b = NA, c = NA, d = NA)),
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

})

