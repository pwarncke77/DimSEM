# Append to tests/testthat/test-bayes-identification.R.
# Tests for facet sign anchoring and CFA-derived starting values.


library(testthat)

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

test_that("dry run reports the weighted constraint; prior reports, never warns", {
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

  # No warning under free + sum_to_zero: the block is genuinely constrained.
  expect_no_warning(
    suppressMessages(.dimsem_estimate_bayes(
      spec, syntax = "", data = dat, covariates = "cov",
      covariate_targets = "all", hyper_loadings = "free",
      covariate_identification = "sum_to_zero",
      bayes_args = list(backend = "dry_run"), verbose = FALSE))
  )

  # Under "prior" the split is left soft BY REQUEST, so the engine reports it
  # rather than warning -- and only when verbose. (Before the weighted
  # projection existed, free loadings raised a warning here; that branch is
  # now unreachable, since every eligible block under "sum_to_zero" is
  # covered by one of the two constraint routes.)
  expect_message(
    fit_p <- .dimsem_estimate_bayes(
      spec, syntax = "", data = dat, covariates = "cov",
      covariate_targets = "all", hyper_loadings = "free",
      covariate_identification = "prior",
      bayes_args = list(backend = "dry_run"), verbose = TRUE),
    "prior-anchored"
  )
  expect_identical(fit_p$stan$data$prj_g, 0L)
  expect_identical(fit_p$stan$identification$weighted, integer(0))
  expect_identical(fit_p$stan$identification$eligible, 1L)

  # verbose = FALSE stays silent on both channels.
  expect_no_warning(
    expect_no_message(.dimsem_estimate_bayes(
      spec, syntax = "", data = dat, covariates = "cov",
      covariate_targets = "all", hyper_loadings = "free",
      covariate_identification = "prior",
      bayes_args = list(backend = "dry_run"), verbose = FALSE))
  )
})


test_that("hyper_scale and hyper_param emit the right Stan under unity", {
  for (gen in list(.dimsem_stan_model, .dimsem_stan_model_gpu)) {
    base <- gen(likelihood = "normal", n_hyper = 1L, n_cov = 1L,
                hyper_loadings = "unity")
    rho  <- gen(likelihood = "normal", n_hyper = 1L, n_cov = 1L,
                hyper_loadings = "unity", hyper_scale = "rho")
    cen  <- gen(likelihood = "normal", n_hyper = 1L, n_cov = 1L,
                hyper_loadings = "unity", hyper_param = "centered")
    both <- gen(likelihood = "normal", n_hyper = 1L, n_cov = 1L,
                hyper_loadings = "unity", hyper_scale = "rho",
                hyper_param = "centered")
    for (code in list(base, rho, cen, both)) {
      ch <- strsplit(code, "")[[1]]
      expect_identical(sum(ch == "{"), sum(ch == "}"))
      # tau_g must be declared exactly once -- as a parameter or as a
      # transformed parameter, never both and never neither.
      n_par <- length(gregexpr("vector<lower=0>\\[G\\] tau_g;", code)[[1]][
        gregexpr("vector<lower=0>\\[G\\] tau_g;", code)[[1]] > 0])
      n_tp  <- length(gregexpr("tau_g = sqrt\\(rho_par", code)[[1]][
        gregexpr("tau_g = sqrt\\(rho_par", code)[[1]] > 0])
      expect_identical(n_par + n_tp, 1L)
    }
    # default is unchanged
    expect_false(grepl("rho_par", base))
    expect_false(grepl("g_resid", base))
    expect_match(base, "tau_g\\[g\\] \\* z_g\\[, g\\]", all = FALSE)

    # rho scale: bounded parameter, derived tau, beta prior, no raw tau prior
    expect_match(rho, "vector<lower=0, upper=1>\\[G\\] rho_par", all = FALSE)
    expect_match(rho, "tau_g = sqrt\\(rho_par", all = FALSE)
    expect_match(rho, "rho_par ~ beta\\(2, 8\\)", all = FALSE)
    expect_false(grepl("tau_g ~ ", rho))

    # centered: disturbances sampled on their own scale
    expect_match(cen, "matrix\\[N, G\\] g_resid", all = FALSE)
    expect_match(cen, "g_resid\\[, g\\] ~ normal\\(0, tau_g\\[g\\]\\)",
                 all = FALSE)
    expect_false(grepl("to_vector\\(z_g\\) ~ std_normal", cen))
    expect_match(cen, "g_score\\[, g\\] = X \\* beta_g\\[, g\\] \\+ g_resid",
                 all = FALSE)

    # both compose
    expect_match(both, "rho_par ~ beta", all = FALSE)
    expect_match(both, "g_resid\\[, g\\] ~ normal\\(0, tau_g\\[g\\]\\)",
                 all = FALSE)

    # custom rho prior is honoured
    r2 <- gen(likelihood = "normal", n_hyper = 1L, n_cov = 1L,
              hyper_loadings = "unity", hyper_scale = "rho",
              priors = list(rho_g = "beta(3, 3)"))
    expect_match(r2, "rho_par ~ beta\\(3, 3\\)", all = FALSE)
  }
})

test_that("the switches are inert under free loadings", {
  for (gen in list(.dimsem_stan_model, .dimsem_stan_model_gpu)) {
    f0 <- gen(likelihood = "normal", n_hyper = 1L, n_cov = 1L,
              hyper_loadings = "free")
    f1 <- gen(likelihood = "normal", n_hyper = 1L, n_cov = 1L,
              hyper_loadings = "free", hyper_scale = "rho",
              hyper_param = "centered")
    expect_identical(f0, f1)          # var(G) = 1: same program either way
    expect_false(grepl("rho_par", f1))
    expect_false(grepl("g_resid", f1))
  }
})

test_that("driver resolves the switches, reports them, and warns under free", {
  set.seed(42)
  n <- 40L
  items <- paste0("it", 1:9)
  dat <- as.data.frame(matrix(rnorm(n * 9), n, 9,
                              dimnames = list(NULL, items)))
  dat$cov <- rnorm(n)
  spec <- list(partition = setNames(rep(1:3, each = 3), items),
               hyper_blocks = list(G1 = c("F1", "F2", "F3")))
  expect_message(
    f <- .dimsem_estimate_bayes(
      spec, syntax = "", data = dat, covariates = "cov",
      covariate_targets = "all", hyper_loadings = "unity",
      covariate_identification = "sum_to_zero",
      bayes_args = list(backend = "dry_run", hyper_scale = "rho"),
      verbose = TRUE),
    "scale = rho"
  )
  expect_match(f$stan$model_code, "rho_par ~ beta", all = FALSE)

  # ignored (with a warning) under free loadings
  expect_warning(
    suppressMessages(.dimsem_estimate_bayes(
      spec, syntax = "", data = dat, covariates = "cov",
      covariate_targets = "all", hyper_loadings = "free",
      covariate_identification = "sum_to_zero",
      bayes_args = list(backend = "dry_run", hyper_param = "centered"),
      verbose = FALSE)),
    "ignored here"
  )
})


test_that("ref_sign recovers relative orientation from the CFA", {
  mk <- function(Psi, latent, cfa_fit = NULL) {
    list(latent_names = latent,
         cfa = list(factor_cor = Psi, fit = cfa_fit))
  }
  # rank-1 correlation structure with known signs; ref items share lavaan's
  # marker orientation (no lavaan object -> o = +1 throughout)
  for (i in 1:20) {
    set.seed(100 + i)
    K <- sample(2:6, 1)
    w <- runif(K, 0.3, 0.9)
    ts <- sample(c(-1, 1), K, replace = TRUE); ts[1] <- 1
    Psi <- outer(w * ts, w * ts); diag(Psi) <- 1
    dimnames(Psi) <- list(paste0("F", 1:K), paste0("F", 1:K))
    got <- suppressWarnings(
      .dimsem_ref_signs(mk(Psi, paste0("F", 1:K)),
                        item_names = paste0("it", seq_len(3 * K)),
                        ref_item = seq(1, 3 * K, by = 3), verbose = FALSE))
    expect_identical(as.numeric(got), as.numeric(ts))
  }
})

test_that("ref_sign degrades safely and warns on an opposed facet", {
  # no CFA at all -> all positive, with a message
  spec <- list(latent_names = paste0("F", 1:3), cfa = NULL)
  expect_message(
    got <- .dimsem_ref_signs(spec, paste0("it", 1:9), c(1L, 4L, 7L),
                             verbose = TRUE),
    "default to positive")
  expect_identical(as.numeric(got), rep(1, 3))

  # genuinely opposed facet -> anchored, but warned about
  Psi <- matrix(c(1, .6, -.5, .6, 1, -.4, -.5, -.4, 1), 3, 3,
                dimnames = list(paste0("F", 1:3), paste0("F", 1:3)))
  expect_warning(
    got <- .dimsem_ref_signs(list(latent_names = paste0("F", 1:3),
                                  cfa = list(factor_cor = Psi, fit = NULL)),
                             paste0("it", 1:9), c(1L, 4L, 7L), verbose = TRUE),
    "oriented OPPOSITE")
  expect_identical(as.numeric(got), c(1, 1, -1))

  # single facet: nothing to orient
  expect_identical(
    as.numeric(.dimsem_ref_signs(list(latent_names = "F1", cfa = NULL),
                                 paste0("it", 1:3), 1L, verbose = FALSE)), 1)
})

test_that("anchored generators reparameterize the loadings correctly", {
  gens <- list(
    function(...) .dimsem_stan_model(likelihood = "normal", ...),
    function(...) .dimsem_stan_model_gpu(likelihood = "normal", ...),
    function(...) .dimsem_stan_model_marginal(...))
  for (gen in gens) {
    off <- gen(n_hyper = 1L, n_cov = 1L, hyper_loadings = "unity",
               anchored = FALSE)
    on  <- gen(n_hyper = 1L, n_cov = 1L, hyper_loadings = "unity",
               anchored = TRUE)
    for (code in list(off, on)) {
      ch <- strsplit(code, "")[[1]]
      expect_identical(sum(ch == "{"), sum(ch == "}"))
      # the anchoring data is declared either way
      expect_match(code, "vector\\[K\\] ref_sign", all = FALSE)
      expect_match(code, "int<lower=0> P_free", all = FALSE)
      expect_match(code, "lam_pos", all = FALSE)
    }
    # off: lambda_un is a plain parameter
    expect_match(off, "vector\\[P\\] lambda_un;", all = FALSE)
    expect_match(off, "lambda_un ~ ", all = FALSE)
    expect_false(grepl("lambda_ref", off))

    # on: split into a constrained reference block and a free remainder
    expect_match(on, "vector<lower=0>\\[K\\] lambda_ref", all = FALSE)
    expect_match(on, "vector\\[P_free\\] lambda_oth", all = FALSE)
    # lambda_un must be DEFINED, not merely referenced: a single
    # declare-and-define built from the anchored pieces via the lam_pos
    # permutation. (The first cut of this shipped the parameters block
    # without the reassembly, leaving lambda_un used but never defined.)
    expect_match(on, "vector\\[P\\] lambda_un =", all = FALSE)
    expect_match(on,
      "append_row\\(ref_sign \\.\\* lambda_ref, lambda_oth\\)\\[lam_pos\\]",
      all = FALSE)
    expect_match(on, "lambda_ref ~ ", all = FALSE)
    expect_match(on, "lambda_oth ~ ", all = FALSE)
    # lambda_un must NOT also be a sampled parameter when anchored
    expect_false(grepl("vector\\[P\\] lambda_un;\\s*//? *(unconstrained|sign fixed)", on))
    # exactly one definition of lambda_un in each mode
    n_def <- function(x) lengths(regmatches(x, gregexpr("vector\\[P\\] lambda_un", x)))
    expect_identical(n_def(off), 1L)
    expect_identical(n_def(on), 1L)
  }
})

test_that("driver defaults anchoring to unity only, and honours the switch", {
  set.seed(42)
  n <- 60L
  items <- paste0("it", 1:9)
  dat <- as.data.frame(matrix(rnorm(n * 9), n, 9,
                              dimnames = list(NULL, items)))
  dat$cov <- rnorm(n)
  spec <- list(partition = setNames(rep(1:3, each = 3), items),
               latent_names = paste0("F", 1:3),
               hyper_blocks = list(G1 = c("F1", "F2", "F3")),
               item_names = items, cfa = NULL)
  run <- function(hl, ...) suppressMessages(suppressWarnings(
    .dimsem_estimate_bayes(spec, syntax = "", data = dat, covariates = "cov",
      covariate_targets = "all", hyper_loadings = hl,
      covariate_identification = "sum_to_zero",
      bayes_args = c(list(backend = "dry_run"), list(...)), verbose = FALSE)))

  expect_match(run("unity")$stan$model_code, "lambda_ref", all = FALSE)   # auto on
  expect_false(grepl("lambda_ref", run("free")$stan$model_code))          # auto off
  expect_false(grepl("lambda_ref", run("unity", anchor_signs = FALSE)$stan$model_code))
  expect_match(run("free", anchor_signs = TRUE)$stan$model_code,
               "lambda_ref", all = FALSE)

  # data is supplied and internally consistent
  d <- run("unity")$stan$data
  expect_identical(length(d$ref_sign), 3L)
  expect_true(all(d$ref_sign %in% c(-1, 1)))
  expect_identical(d$P_free, 9L - 3L)
  expect_identical(sort(c(d$ref_item, d$nonref_item)), 1:9)
  # lam_pos is a genuine permutation, and it puts the reference items first
  expect_identical(sort(d$lam_pos), 1:9)
  expect_identical(d$lam_pos[d$ref_item], 1:3)
  expect_true(all(d$lam_pos[d$nonref_item] > 3L))
})
