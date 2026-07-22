test_that("method normalization expands all", {
  expect_equal(
    sort(DimSem:::.dimsem_normalize_methods("all")),
    c("EGA", "PA", "TNN")
  )
})
