source("blavaan_comparison.R")
ok <- function(l) cat("PASS:", l, "\n")

fac <- paste0("Fat_", 1:6)
raw <- c(paste0("ideology=~", fac), paste0(fac, "~lib_con"), "ideology~lib_con",
         paste0(fac, "~~", fac), "ideology~~ideology")

# 1. exact
for (f in fac) stopifnot(raw[.dimsem_match_par(raw, "ideology", "=~", f)] ==
                           paste0("ideology=~", f))
stopifnot(raw[.dimsem_match_par(raw, "ideology", "~", "lib_con")] == "ideology~lib_con")
for (f in fac) stopifnot(raw[.dimsem_match_par(raw, f, "~", "lib_con")] ==
                           paste0(f, "~lib_con"))
ok("exact match; '~' vs '=~' vs '~~' never confused")

# 2. whitespace variants
sp <- gsub("~", " ~ ", raw); sp <- gsub("= ~ ", "=~", sp)
for (f in fac) stopifnot(!is.na(.dimsem_match_par(sp, f, "~", "lib_con")))
ok("whitespace-insensitive fallback")

# 3. make.names() mangling (read.csv round trip)
mang <- make.names(raw, unique = TRUE)
for (f in fac) {
  i <- .dimsem_match_par(mang, "ideology", "=~", f)
  stopifnot(mang[i] == make.names(paste0("ideology=~", f)))
}
ok("punctuation-insensitive fallback recovers make.names() mangling")

# 4. failure modes are informative, not silent
e <- tryCatch(.dimsem_match_par(raw, "ideology", "~", "income"),
              error = function(e) conditionMessage(e))
stopifnot(is.character(e), grepl("Could not locate", e), grepl("ideology", e))
ok("missing parameter -> informative error listing near matches")

# 5. end-to-end on a mock standardizedPosterior
set.seed(1); n <- 500
d <- matrix(rnorm(n * length(raw)), n, length(raw), dimnames = list(NULL, raw))
delta_true <- c(.80, .64, .83, .72, .52, .45)
alpha <- rnorm(n, .5, .12)
for (k in seq_along(fac)) {
  d[, paste0("ideology=~", fac[k])] <- delta_true[k] + rnorm(n, 0, .02)
  d[, paste0(fac[k], "~lib_con")] <- rnorm(n, .05, .1) - delta_true[k] * alpha
}
d[, "ideology~lib_con"] <- alpha
local({
  syntax <- list(hyper_name = "ideology", factors = fac, covariate = "lib_con")
  # Shim the `::` operator so blavaan::standardizedPosterior resolves to the
  # mock draws without the package being installed.
  shim <- new.env(parent = globalenv())
  assign("::", function(pkg, name) {
    pkg <- as.character(substitute(pkg)); name <- as.character(substitute(name))
    if (pkg == "blavaan" && name == "standardizedPosterior") return(function(fit) d)
    get(name, envir = asNamespace(pkg))
  }, envir = shim)
  assign("requireNamespace", function(...) TRUE, envir = shim)

  f <- blavaan_covariate_decomposition; environment(f) <- shim
  out <- f(NULL, syntax)
  stopifnot(nrow(out) == 1 + 6 + 6 + 1 + 6 + 6)
  b <- out[grepl("sum-to-zero", out$parameter, fixed = TRUE) &
             startsWith(out$parameter, "Fat_"), ]
  cat("  projected betas:", sprintf("%+.3f", b$est), "\n")
  rd <- attr(out, "ridge_diagnostics")
  stopifnot(all(rd$cor_alpha_beta < -0.4), all(rd$ratio_sd > 1))
  cat("  corr(alpha, beta):", sprintf("%.2f", rd$cor_alpha_beta), "\n")
  tot <- out[grepl("(total)", out$parameter, fixed = TRUE), ]
  cat("  sd(total) vs sd(beta):", sprintf("%.3f", tot$se), "|",
      sprintf("%.3f", rd$sd_beta), "\n")

  g <- blavaan_hyper_strength; environment(g) <- shim
  hs <- g(NULL, syntax)
  truth <- apply(combn(6, 2), 2, function(p) delta_true[p[1]] * delta_true[p[2]])
  stopifnot(nrow(hs) == 15, max(abs(hs$implied_r - truth)) < 0.01)
  cat("  implied pair r range:", sprintf("%.3f-%.3f", min(hs$implied_r),
      max(hs$implied_r)), "\n")
})
ok("end-to-end: decomposition, ridge diagnostics, implied pair correlations")
cat("\nAll matcher checks passed.\n")
