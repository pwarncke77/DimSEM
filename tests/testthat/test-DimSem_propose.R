test_that("DimSem_propose returns a basic PA-only proposal", {
  skip_if_not_installed("psych")
  skip_if_not_installed("igraph")
  skip_if_not_installed("ggraph")
  skip_if_not_installed("tidygraph")
  skip_if_not_installed("ggplot2")

  set.seed(1)
  x <- matrix(rnorm(150), nrow = 30, ncol = 5)
  colnames(x) <- paste0("i", 1:5)

  out <- DimSem_propose(
    x,
    methods = "PA",
    partition_source = "PA",
    pa_args = list(n.iter = 5, plot = FALSE),
    make_plots = FALSE,
    verbose = FALSE
  )

  expect_s3_class(out, "DimSem_proposal")
  expect_true(is.data.frame(out$selected$partition_table))
  expect_equal(nrow(out$selected$partition_table), 5)
  expect_true(is.character(out$selected$model_syntax))
})
