context("parse_delimitted_string tests")
library(miscutil)

test_that("parse delimitter", {
  x = "one,two,three,"
  expect_equal(length(parse_delimitted_string(x)),3)
  expect_equal(length(parse_delimitted_string(x, rm.empty = F)),4)
  expect_equal(length(parse_delimitted_string(x, delim = 'w', rm.empty = F)),2)
  expect_match(parse_delimitted_string(x)[1],'one')
  expect_is(parse_delimitted_string(x), 'character')
})
