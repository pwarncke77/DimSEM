
test_that("free loadings + sum_to_zero take the weighted-projection route", {
  bb <- .dimsem_beta_basis(c(1L, 1L, 1L, 1L), factor_hyper = c(1L, 1L, 1L, 1L),
                           reg_g = 1L, hyper_loadings = "free",
                           covariate_identification = "sum_to_zero")
  expect_identical(bb$eligible, 1L)
  expect_identical(bb$constrained, integer(0))   # no fixed basis under free
  expect_identical(bb$weighted, 1L)              # runtime projection instead
  expect_identical(bb$K_free, 4L)                # identity coords, m of them
  expect_equal(bb$basis, diag(4), ignore_attr = TRUE)

  # unity keeps the fixed Helmert basis, and never the projection
  bu <- .dimsem_beta_basis(c(1L, 1L, 1L, 1L), factor_hyper = c(1L, 1L, 1L, 1L),
                           reg_g = 1L, hyper_loadings = "unity",
                           covariate_identification = "sum_to_zero")
  expect_identical(bu$constrained, 1L)
  expect_identical(bu$weighted, integer(0))
  expect_identical(bu$K_free, 3L)

  # "prior" takes neither route
  bp <- .dimsem_beta_basis(c(1L, 1L, 1L, 1L), factor_hyper = c(1L, 1L, 1L, 1L),
                           reg_g = 1L, hyper_loadings = "free",
                           covariate_identification = "prior")
  expect_identical(bp$constrained, integer(0))
  expect_identical(bp$weighted, integer(0))
})

test_that("generators emit the projection under free + sum_to_zero only", {
  for (gen in list(.dimsem_stan_model, .dimsem_stan_model_gpu)) {
    for (lik in c("normal", "ordinal")) {
      fr <- gen(likelihood = lik, n_hyper = 1L, n_cov = 1L,
                hyper_loadings = "free")
      un <- gen(likelihood = lik, n_hyper = 1L, n_cov = 1L,
                hyper_loadings = "unity")
      for (code in list(fr, un)) {
        ch <- strsplit(code, "")[[1]]
        expect_identical(sum(ch == "{"), sum(ch == "}"))
        expect_match(code, "array\\[G\\] int<lower=0, upper=1> prj_g",
                     all = FALSE)
      }
      # the projection body is free-only; the fixed-basis relabeling unity-only
      expect_match(fr, "prj_g\\[g\\] == 1", all = FALSE)
      expect_false(grepl("prj_g\\[g\\] == 1", un))
      expect_match(fr, "beta_f\\[q, k\\] -= delta_l\\[k\\] \\* num / den",
                   all = FALSE)
      expect_match(un, "ctr_g\\[g\\] == 1", all = FALSE)
      # projection must come after delta_l is populated and before eta
      pos_d <- regexpr("delta_l\\[k\\] = factor_hyper", fr)
      pos_p <- regexpr("prj_g\\[g\\] == 1", fr)
      pos_e <- regexpr("eta\\[, k\\] =", fr)
      expect_lt(pos_d, pos_p)
      expect_lt(pos_p, pos_e)
      # no covariates -> no projection machinery
      expect_false(grepl("prj_g\\[g\\] == 1",
                         gen(likelihood = lik, n_hyper = 1L, n_cov = 0L,
                             hyper_loadings = "free")))
    }
  }
})

test_that("weighted projection satisfies its constraint and reduces to unity", {
  # R replica of the Stan projection.
  project <- function(beta, delta) {
    den <- sum(delta^2)
    if (den <= 1e-12) return(beta)
    beta - delta * sum(delta * beta) / den
  }
  set.seed(11)
  d <- c(2.834, 0.869, 0.473, 0.200)
  for (i in 1:200) {
    g <- rnorm(4, 0, 0.5)
    b <- project(g, d)
    expect_lt(abs(sum(d * b)), 1e-12)          # constraint exact
    alpha <- rnorm(1)
    t <- b + d * alpha
    expect_equal(sum(d * t) / sum(d^2), alpha, tolerance = 1e-12)
  }
  # sign-orbit invariance: flipping (delta_k, beta_k) jointly
  for (bits in 0:15) {
    s <- ifelse(bitwAnd(bits, 2^(0:3)) > 0, -1, 1)
    g <- rnorm(4, 0, 0.5)
    expect_equal(project(g * s, d * s), project(g, d) * s, tolerance = 1e-12)
  }
  # delta == 1 reduces to mean-centring
  g <- rnorm(4)
  expect_equal(project(g, rep(1, 4)), g - mean(g), tolerance = 1e-12)
  # degenerate block: zero loadings leave the coefficients untouched
  expect_equal(project(g, rep(0, 4)), g, tolerance = 1e-12)
})

test_that("dry run reports the weighted constraint without the soft warning", {
  set.seed(42)
  n <- 40L
  items <- paste0("it", 1:9)
  dat <- as.data.frame(matrix(rnorm(n * 9), n, 9,
                              dimnames = list(NULL, items)))
  dat$cov <- rnorm(n)
  spec <- list(partition = setNames(rep(1:3, each = 3), items),
               hyper_blocks = list(G1 = c("F1", "F2", "F3")))

  expect_message(
    fit <- .dimsem_estimate_bayes(
      spec, syntax = "", data = dat, covariates = "cov",
      covariate_targets = "all", hyper_loadings = "free",
      covariate_identification = "sum_to_zero",
      bayes_args = list(backend = "dry_run"), verbose = TRUE),
    "delta-weighted sum-to-zero"
  )
  expect_identical(fit$stan$data$prj_g, 1L)
  expect_identical(fit$stan$data$ctr_g, 0L)
  expect_identical(fit$stan$identification$weighted, 1L)

  # the free-loadings soft-identification warning must NOT fire any more
  expect_no_warning(
    suppressMessages(.dimsem_estimate_bayes(
      spec, syntax = "", data = dat, covariates = "cov",
      covariate_targets = "all", hyper_loadings = "free",
      covariate_identification = "sum_to_zero",
      bayes_args = list(backend = "dry_run"), verbose = FALSE))
  )
  # but it still fires under "prior"
  expect_warning(
    suppressMessages(.dimsem_estimate_bayes(
      spec, syntax = "", data = dat, covariates = "cov",
      covariate_targets = "all", hyper_loadings = "free",
      covariate_identification = "prior",
      bayes_args = list(backend = "dry_run"), verbose = FALSE)),
    "identified only softly"
  )
})
