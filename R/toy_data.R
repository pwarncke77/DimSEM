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

# --- write to data/ ------------------------------------------------------------

# Equivalent to usethis::use_data(..., overwrite = TRUE):
save(toy_simple,  file = "data/toy_simple.rda",  compress = "bzip2", version = 2)
save(toy_complex, file = "data/toy_complex.rda", compress = "bzip2", version = 2)
save(toy_random,  file = "data/toy_random.rda",  compress = "bzip2", version = 2)

