# =============================================================================
# data-raw/toy_data.R
# -----------------------------------------------------------------------------
# Generates the three toy datasets shipped with DimSem:
#
#   toy_simple  : n = 400, 12 five-point Likert items, 3 latent factors
#                 (r = .30 between all factor pairs). Items a:d -> F1,
#                 e:h -> F2, i:l -> F3.
#
#   toy_complex : n = 400, 24 five-point Likert items. Items a:j -> F1
#                 (10 items), k:t -> F2 (10 items), u:w -> F3 (3 items),
#                 item x is an orphan (pure noise, unrelated to anything).
#                 cor(F1, F2) = .40; F3 uncorrelated with F1 and F2.
#
#   toy_random  : n = 400, 24 five-point Likert items, no factor structure
#                 (independent items).
#
#   toy_simple_cov : n = 400, toy_simple's 3-factor/12-item structure plus a
#                 latent-level external covariate "cov" with known ground
#                 truth for the sum-to-zero covariate decomposition
#                 (see section 4; truth in attr(toy_simple_cov, "truth")).
#
# Run this script from the package root. It writes .rda files into data/.
# The script itself lives in data-raw/ and is excluded from the build
# (add ^data-raw$ to .Rbuildignore).
#
# Base R only -- no package dependencies.
# =============================================================================

set.seed(123)

# --- helpers -----------------------------------------------------------------

# Draw n multivariate-normal rows with mean 0 and covariance Sigma.
rmvn <- function(n, Sigma) {
  p <- ncol(Sigma)
  matrix(stats::rnorm(n * p), nrow = n, ncol = p) %*% chol(Sigma)
}

# Discretize a standardized continuous variable into a 5-point Likert item.
# Symmetric thresholds yield an approximately bell-shaped category
# distribution centered on the scale midpoint (3).
to_likert5 <- function(x) {
  thresholds <- c(-1.5, -0.5, 0.5, 1.5)
  as.integer(findInterval(x, thresholds) + 1L)
}

# Generate a Likert-type item set from a common-factor model.
#   n         : sample size
#   Phi       : factor correlation matrix (k x k)
#   assign    : integer vector of length p mapping each item to a factor
#               (NA = orphan item, i.e., pure noise)
#   loadings  : numeric vector of length p (standardized loadings; ignored
#               for orphan items)
#   item_names: character vector of length p
make_factor_data <- function(n, Phi, assign, loadings, item_names) {
  p <- length(assign)
  stopifnot(length(loadings) == p, length(item_names) == p)

  FF <- rmvn(n, Phi)                          # latent factor scores
  X  <- matrix(NA_real_, nrow = n, ncol = p)

  for (j in seq_len(p)) {
    if (is.na(assign[j])) {
      # Orphan item: pure standard-normal noise.
      X[, j] <- stats::rnorm(n)
    } else {
      lam    <- loadings[j]
      X[, j] <- lam * FF[, assign[j]] + sqrt(1 - lam^2) * stats::rnorm(n)
    }
  }

  out <- as.data.frame(apply(X, 2, to_likert5))
  names(out) <- item_names
  out
}

# --- 1) toy_simple: 3 correlated factors, 4 items each ------------------------

Phi_simple <- matrix(0.4, nrow = 3, ncol = 3)
diag(Phi_simple) <- 1

toy_simple <- make_factor_data(
  n          = 400,
  Phi        = Phi_simple,
  assign     = rep(1:3, each = 4),                 # a:d, e:h, i:l
  loadings   = rep(c(0.80, 0.75, 0.70, 0.65), 3),  # mixed loading strengths
  item_names = letters[1:12]
)

# --- 2) toy_complex: 10 + 10 + 3 items plus one orphan ------------------------

Phi_complex <- diag(3)
Phi_complex[1, 2] <- Phi_complex[2, 1] <- 0.4      # F1-F2 = .40; F3 free-standing

toy_complex <- make_factor_data(
  n          = 400,
  Phi        = Phi_complex,
  assign     = c(rep(1L, 10), rep(2L, 10), rep(3L, 3), NA),
  loadings   = c(
    rep(c(0.80, 0.75, 0.70, 0.65, 0.60), 2),       # F1: items a:j
    rep(c(0.80, 0.75, 0.70, 0.65, 0.60), 2),       # F2: items k:t
    c(0.80, 0.75, 0.70),                           # F3: items u:w
    NA                                             # orphan: item x
  ),
  item_names = letters[1:24]
)

# --- 3) toy_random: no factor structure ---------------------------------------

toy_random <- as.data.frame(
  apply(matrix(stats::rnorm(400 * 24), nrow = 400), 2, to_likert5)
)
names(toy_random) <- letters[1:24]

# --- 4) toy_simple_cov: toy_simple's structure plus a latent-level covariate --
#
# Validation target for simultaneous covariate estimation under unity
# hyper-loadings with the sum-to-zero identification. The covariate is
# generated at the LATENT level (never from items), so the measurement
# model is correctly specified: items are conditionally independent given
# (eta, cov), and the residual (conditional-on-cov) factor correlation is
# exactly compound-symmetric -- the healthy-tau regime in which the
# sign-relabeling block premise holds.
#
# Generating model, in the estimator's own metric (disturbance var = 1):
#   cov     ~ N(0, 1)
#   g_resid ~ N(0, tau2),      tau2 = 0.43  =>  rho_g = tau2/(tau2+1) = 0.30
#   eta_k   = t_k * cov + g_resid + z_k,   z_k ~ N(0, 1) iid
# with total covariate effects t = (0.60, 0.30, 0.45). Under the
# sum-to-zero convention the true decomposition is therefore
#   alpha = mean(t) = 0.45,   beta = t - alpha = (0.15, -0.15, 0.00).
# Items are built from the standardized eta with the same loading pattern
# as toy_simple and discretized to 5-point Likert. Ground truth is exact
# at the latent level; Likert coarsening mildly attenuates estimates
# under the normal likelihood, so validation should target coverage and
# ordering of the recovered decomposition, not equality.
#
# The truth is attached as attr(toy_simple_cov, "truth"). Note that
# DimSEM_propose() should be run on the item columns only, e.g.
# DimSEM_propose(toy_simple_cov[, letters[1:12]]).

set.seed(456)

n_cov      <- 400
t_totals   <- c(0.60, 0.30, 0.45)
tau2_cov   <- 0.43
cov_scores <- stats::rnorm(n_cov)
g_resid    <- stats::rnorm(n_cov, sd = sqrt(tau2_cov))
eta_cov    <- sapply(1:3, function(k) {
  t_totals[k] * cov_scores + g_resid + stats::rnorm(n_cov)
})

loadings_cov <- rep(c(0.80, 0.75, 0.70, 0.65), 3)
assign_cov   <- rep(1:3, each = 4)
X_cov <- sapply(seq_len(12), function(j) {
  lam <- loadings_cov[j]
  eta_std <- as.numeric(scale(eta_cov[, assign_cov[j]]))
  lam * eta_std + sqrt(1 - lam^2) * stats::rnorm(n_cov)
})

toy_simple_cov <- as.data.frame(apply(X_cov, 2, to_likert5))
names(toy_simple_cov) <- letters[1:12]
toy_simple_cov$cov <- cov_scores
attr(toy_simple_cov, "truth") <- list(
  totals = t_totals,
  alpha  = mean(t_totals),
  beta   = t_totals - mean(t_totals),
  tau2   = tau2_cov,
  rho_g  = tau2_cov / (tau2_cov + 1)
)

# Equivalent to usethis::use_data(..., overwrite = TRUE):
save(toy_simple,  file = "data/toy_simple.rda",  compress = "bzip2", version = 2)
save(toy_complex, file = "data/toy_complex.rda", compress = "bzip2", version = 2)
save(toy_random,  file = "data/toy_random.rda",  compress = "bzip2", version = 2)
save(toy_simple_cov, file = "data/toy_simple_cov.rda",
     compress = "bzip2", version = 2)

