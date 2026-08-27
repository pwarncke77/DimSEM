# Append to tests/testthat/test-bayes-identification.R (or its own file).
# Tests for the marginalized Stan program.

test_that("marginal generator emits a complete, well-formed program", {
  for (hl in c("unity", "free")) {
    code <- .dimsem_stan_model_marginal(n_hyper = 1L, n_cov = 1L,
                                        hyper_loadings = hl)
    ch <- strsplit(code, "")[[1]]
    expect_identical(sum(ch == "{"), sum(ch == "}"))
    for (b in c("data \\{", "transformed data \\{", "parameters \\{",
                "transformed parameters \\{", "model \\{",
                "generated quantities \\{")) {
      expect_match(code, b, all = FALSE)
    }
    # pattern machinery
    for (v in c("pat_r_start", "pat_r_end", "pat_i_start", "pat_i_len",
                "pat_items", "n_pat_items")) {
      expect_match(code, v, all = FALSE)
    }
    # Woodbury pieces
    expect_match(code, "cholesky_decompose", all = FALSE)
    expect_match(code, "mdivide_left_tri_low", all = FALSE)
    expect_match(code, "rows_dot_self", all = FALSE)
    expect_match(code, "columns_dot_self", all = FALSE)
    expect_match(code, "int R = K \\+ G", all = FALSE)
    # the conditional program's long format must be gone
    expect_false(grepl("n_obs", code))
    expect_false(grepl("reduce_sum", code))
    # reporting survives unchanged
    for (v in c("beta_f_cor", "beta_g_cor", "beta_total", "lambda", "sign_f")) {
      expect_match(code, v, all = FALSE)
    }
    # loading-mode specifics
    if (hl == "unity") {
      expect_match(code, "lambda_un\\[p\\] \\* tau_g\\[g\\]", all = FALSE)
      expect_match(code, "rho_g", all = FALSE)
    } else {
      expect_match(code, "lambda_un\\[p\\] \\* delta_l\\[item_factor\\[p\\]\\]",
                   all = FALSE)
      expect_match(code, "prj_g\\[g\\] == 1", all = FALSE)
    }
  }
  # rho scale composes
  r <- .dimsem_stan_model_marginal(n_hyper = 1L, n_cov = 1L,
                                   hyper_loadings = "unity",
                                   hyper_scale = "rho")
  expect_match(r, "vector<lower=0, upper=1>\\[G\\] rho_par", all = FALSE)
  expect_match(r, "tau_g = sqrt\\(rho_par", all = FALSE)
  # no covariates / no hyper still produce valid programs
  for (nh in 0:1) for (nc in 0:1) {
    cc <- .dimsem_stan_model_marginal(n_hyper = nh, n_cov = nc)
    ch <- strsplit(cc, "")[[1]]
    expect_identical(sum(ch == "{"), sum(ch == "}"))
  }
})

test_that("pattern builder sorts respondents and indexes items correctly", {
  set.seed(1)
  N <- 30L; P <- 5L
  Y <- matrix(rnorm(N * P), N, P)
  Y[1:5, 2] <- NA          # pattern A
  Y[6:10, c(2, 4)] <- NA   # pattern B
  Y[11, ] <- NA            # no observed items -> dropped
  sd0 <- list(N = N, P = P, K = 2L, Y_raw = Y,
              item_factor = c(1L, 1L, 2L, 2L, 2L),
              n_obs = 1L, obs_row = 1L, obs_col = 1L, y_obs = 1,
              grainsize = 1L, X = matrix(rnorm(N), N, 1L))
  out <- suppressMessages(.dimsem_stan_data_marginal(sd0))

  expect_identical(out$N, N - 1L)                    # the all-NA row is gone
  expect_null(out$n_obs); expect_null(out$Y_raw)     # long format removed
  expect_identical(nrow(out$X), out$N)
  expect_true(all(is.finite(out$Y)))                 # NAs zero-filled

  # every pattern's row range is contiguous and the observed items match
  for (s in seq_len(out$n_pat)) {
    idx <- out$pat_items[out$pat_i_start[s] +
                         seq_len(out$pat_i_len[s]) - 1L]
    rows <- out$pat_r_start[s]:out$pat_r_end[s]
    expect_gt(length(idx), 0L)
    # all rows in a pattern share exactly this observed-item set
    for (r in rows) {
      obs_r <- which(out$Y[r, ] != 0 | seq_len(P) %in% idx)
      expect_true(all(idx %in% obs_r))
    }
  }
  expect_identical(sum(out$pat_i_len), out$n_pat_items)
  expect_identical(sum(out$pat_r_end - out$pat_r_start + 1L), out$N)
})

test_that("marginal switch is refused where it cannot apply", {
  set.seed(42)
  n <- 40L
  items <- paste0("it", 1:9)
  dat <- as.data.frame(matrix(rnorm(n * 9), n, 9,
                              dimnames = list(NULL, items)))
  dat$cov <- rnorm(n)
  spec <- list(partition = setNames(rep(1:3, each = 3), items),
               hyper_blocks = list(G1 = c("F1", "F2", "F3")))
  base <- list(spec, syntax = "", data = dat, covariates = "cov",
               covariate_targets = "all")

  expect_error(
    .dimsem_estimate_bayes(spec, syntax = "", data = dat, covariates = "cov",
      covariate_targets = "all", hyper_loadings = "unity",
      bayes_args = list(backend = "dry_run", marginal = TRUE,
                        likelihood = "ordinal"), verbose = FALSE),
    "requires likelihood"
  )
  expect_error(
    .dimsem_estimate_bayes(spec, syntax = "", data = dat, covariates = "cov",
      covariate_targets = "all", hyper_loadings = "unity",
      bayes_args = list(backend = "dry_run", marginal = TRUE, gpu = TRUE),
      verbose = FALSE),
    "mutually exclusive"
  )
  # and it works where it does apply
  f <- suppressMessages(.dimsem_estimate_bayes(
    spec, syntax = "", data = dat, covariates = "cov",
    covariate_targets = "all", hyper_loadings = "unity",
    covariate_identification = "sum_to_zero",
    bayes_args = list(backend = "dry_run", marginal = TRUE), verbose = FALSE))
  expect_match(f$stan$model_code, "pat_r_start", all = FALSE)
  expect_true(!is.null(f$stan$data$n_pat))
  expect_identical(f$stan$data$ctr_g, 1L)   # identification still wired
})
