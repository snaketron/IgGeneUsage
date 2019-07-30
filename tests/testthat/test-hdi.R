context("Tests HDI")


test_that("HDI calculation check", {
  
  expect_error(object = getHdi(vec = "a", hdi.level = 0.95),
               regexp = "HDI not numeric")
  
  expect_is(object = getHdi(vec = 1, hdi.level = 0.95),
            class = "numeric")
  
  hdi <- getHdi(vec = rnorm(n = 10^6, mean = 0, sd = 1), hdi.level = 0.95)
  expect_gte(object = hdi[1], expected = -2)
  expect_lte(object = hdi[2], expected = 2)
  
  hdi <- getHdi(vec = rnorm(n = 10^6, mean = 0, sd = 100), hdi.level = 0.95)
  expect_lte(object = hdi[1], expected = -190)
  expect_gte(object = hdi[2], expected = 190)
})
