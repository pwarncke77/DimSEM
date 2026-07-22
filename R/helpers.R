# Internal utilities for dimsem -------------------------------------------------
#
# This file is the single canonical home for all .dimsem_* internal helpers.
# Do not redefine any of these functions in other R/ files: files in R/ are
# sourced into one shared namespace in collation order, so a duplicate
# definition in a later file silently overrides an earlier one.

# --- data preparation ----------------------------------------------------------

# Validate `data`, optionally subset item columns via `items`, and optionally
# rename the analyzed items via `item_names`.
#
#   items      : NULL (use all columns), a character vector of column names,
#                or an integer vector of column positions. Subsetting happens
#                *before* the numeric check, so non-numeric ID columns outside
#                the selection are allowed.
#   item_names : NULL (use column names, or item_1, item_2, ... when absent),
#                or a character vector with one name per *analyzed* item.
.dimsem_prepare_data <- function(data, items = NULL, item_names = NULL) {
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("`data` must be a data frame or matrix.", call. = FALSE)
  }

  data <- as.data.frame(data)

  if (!is.null(items)) {
    if (is.character(items)) {
      if (is.null(colnames(data))) {
        stop("`items` is a character vector, but `data` has no column names.",
             call. = FALSE)
      }
      missing_items <- setdiff(items, colnames(data))
      if (length(missing_items) > 0) {
        stop("`items` not found in `data`: ",
             paste(missing_items, collapse = ", "), call. = FALSE)
      }
    } else if (is.numeric(items)) {
      items <- as.integer(items)
      if (anyNA(items) || any(items < 1L) || any(items > ncol(data))) {
        stop("Integer `items` must be column positions between 1 and ",
             ncol(data), ".", call. = FALSE)
      }
    } else {
      stop("`items` must be a character vector of column names or an ",
           "integer vector of column positions.", call. = FALSE)
    }

    if (anyDuplicated(items)) {
      stop("`items` must not contain duplicated entries.", call. = FALSE)
    }
    data <- data[, items, drop = FALSE]
  }

  numeric_cols <- vapply(data, is.numeric, logical(1))
  if (!all(numeric_cols)) {
    bad <- paste(names(data)[!numeric_cols], collapse = ", ")
    stop("All analyzed item columns must be numeric. Non-numeric columns: ",
         bad, call. = FALSE)
  }

  p <- ncol(data)
  if (p < 2) {
    stop("`data` must contain at least two item columns after applying ",
         "`items`.", call. = FALSE)
  }

  if (is.null(item_names)) {
    item_names <- colnames(data)
    if (is.null(item_names) || anyNA(item_names) || any(item_names == "")) {
      item_names <- paste0("item_", seq_len(p))
    }
  } else {
    item_names <- as.character(item_names)
  }

  if (length(item_names) != p) {
    stop("`item_names` must have length equal to the number of analyzed ",
         "items (", p, "), but has length ", length(item_names), ".",
         call. = FALSE)
  }
  if (anyNA(item_names) || any(item_names == "")) {
    stop("`item_names` must not contain missing or empty names.",
         call. = FALSE)
  }
  if (anyDuplicated(item_names)) {
    stop("`item_names` must be unique.", call. = FALSE)
  }

  colnames(data) <- item_names
  list(data = data, item_names = item_names)
}

.dimsem_normalize_methods <- function(methods) {
  if (length(methods) == 0L) {
    stop("`methods` must contain at least one method.", call. = FALSE)
  }

  methods <- trimws(toupper(methods))
  if ("ALL" %in% methods) {
    methods <- c("EGA", "PA", "EKC")
  }

  allowed <- c("EGA", "PA", "EKC")
  bad <- setdiff(methods, allowed)
  if (length(bad) > 0) {
    stop("Unknown method(s): ", paste(bad, collapse = ", "), call. = FALSE)
  }

  unique(methods)
}

# --- method runners -------------------------------------------------------------

# Level-1 EGA with a minimum-community-size rule (default >= 3 items).
#
# Rationale: a two-item community cannot identify a latent dimension on its
# own (the latent-variable analogue of the three-indicator rule; Gorsuch,
# 1997), mirroring both the two-node-community handling in Taxonomic Graph
# Analysis (Samo et al., 2025) and dimsem's own hyper-level
# `min_block_size` rule. Enforcement uses resolution annealing: EGA is run
# with the Louvain algorithm, and whenever the solution contains fixable
# sub-threshold communities, the Louvain resolution parameter is lowered
# stepwise (lower resolution = coarser partition) until every community has
# at least `min_community_size` items -- or the resolution floor is
# reached, in which case remaining sub-threshold communities are converted
# to unassigned orphans (NA) with a warning.
#
# A sub-threshold community is "fixable" by coarsening only if it has at
# least one edge to the rest of the network. Communities that form their
# own disconnected component (an isolated item, or an isolated doublet)
# can NEVER be merged by lowering the resolution -- Louvain cannot join
# disconnected components -- so they are exempted from the annealing loop
# and orphaned directly, preventing a futile walk to the resolution floor.
#
# Intercepted convenience arguments in `ega_args` (not passed to
# EGAnet::EGA): `min_community_size` (default 3; set to 1 to disable the
# rule), `resolution_start` (default 1.2), `resolution_step` (default
# 0.01), `resolution_min` (default 0.05). Supplying `resolution` directly
# in `ega_args` fixes the resolution and disables the search (sub-threshold
# communities are then orphaned immediately). The rule requires
# `algorithm = "louvain"` (the new default, with
# `consensus.method = "lowest_tefi"`); overriding to another algorithm
# (e.g., "walktrap") also disables the search.
.dimsem_run_ega <- function(data, ega_args = list(), seed = NULL,
                            verbose = TRUE) {
  .dimsem_require_namespace("EGAnet")

  defaults <- list(
    corr = "auto",
    model = "glasso",
    algorithm = "louvain",
    consensus.method = "lowest_tefi",
    uni.method = "LE",
    plot.EGA = FALSE
  )
  args <- .dimsem_merge_args(defaults, ega_args)

  # Intercept convenience arguments.
  min_size <- args$min_community_size %||% 3L
  res_start <- args$resolution_start %||% 1.2
  res_step <- args$resolution_step %||% 0.01
  res_min <- args$resolution_min %||% 0.05
  fixed_res <- args$resolution
  args$min_community_size <- NULL
  args$resolution_start <- NULL
  args$resolution_step <- NULL
  args$resolution_min <- NULL
  args$resolution <- NULL

  louvain <- identical(tolower(as.character(args$algorithm %||% "")),
                       "louvain")
  search_enabled <- louvain && min_size > 1L && is.null(fixed_res)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  run_once <- function(resolution = NULL) {
    call_args <- args
    if (louvain && !is.null(resolution)) {
      call_args$resolution <- resolution
    }

    ega_empty_warning <- FALSE
    withCallingHandlers(
      do.call(EGAnet::EGA, c(list(data = data), call_args)),
      warning = function(w) {
        msg <- conditionMessage(w)

        if (grepl("network.*input.*empty|input.*network.*empty|network.*empty",
                  msg, ignore.case = TRUE)) {
          ega_empty_warning <<- TRUE
          invokeRestart("muffleWarning")
        }

        if (isTRUE(ega_empty_warning) &&
            grepl("Some variables did not belong|TEFI calculation",
                  msg, ignore.case = TRUE)) {
          invokeRestart("muffleWarning")
        }
      }
    )
  }

  resolution <- if (louvain) fixed_res %||% res_start else NULL
  n_steps <- 0L
  # Best-solution tracking: annealing does not improve monotonically --
  # near the resolution floor Louvain can collapse the entire network into
  # one community, destroying genuine structure while an unmergeable
  # doublet persists. The accepted solution is therefore the one with the
  # FEWEST sub-threshold items along the whole trajectory, with ties
  # broken toward the HIGHEST resolution (finest structure-preserving
  # partition); the search itself also stops once a solution coarsens to a
  # single community, since nothing coarser exists.
  best <- NULL

  repeat {
    ega <- run_once(resolution)
    wc <- .dimsem_extract_ega_wc(ega)
    names(wc) <- colnames(data)
    network <- .dimsem_extract_ega_network(ega, item_names = colnames(data))

    sizes <- table(wc[!is.na(wc)])
    small <- names(sizes)[sizes < min_size]
    n_small_items <- if (length(small) > 0) {
      sum(wc %in% as.integer(small), na.rm = TRUE)
    } else {
      0L
    }

    if (is.null(best) || n_small_items < best$n_small_items) {
      best <- list(ega = ega, wc = wc, network = network,
                   resolution = resolution,
                   n_small_items = n_small_items)
    }

    # Coarsening can only help sub-threshold communities that have some
    # connection to the rest of the network; communities forming their own
    # disconnected component can never be merged by lowering resolution.
    fixable <- FALSE
    if (length(small) > 0 && !is.null(network)) {
      W <- abs(network)
      diag(W) <- 0
      fixable <- any(vapply(small, function(cc) {
        members <- names(wc)[!is.na(wc) & wc == as.integer(cc)]
        others <- setdiff(colnames(W), members)
        length(others) > 0 &&
          sum(W[members, others, drop = FALSE]) > .Machine$double.eps
      }, logical(1)))
    }

    single_community <- length(sizes) <= 1
    can_step <- search_enabled && fixable && !single_community &&
      n_small_items > 0L &&
      !is.null(resolution) && (resolution - res_step) >= res_min
    if (!can_step) {
      break
    }
    resolution <- resolution - res_step
    n_steps <- n_steps + 1L
  }

  # Accept the best solution seen (identical to the last one whenever the
  # search ended cleanly with zero sub-threshold items).
  ega <- best$ega
  wc <- best$wc
  network <- best$network
  accepted_resolution <- best$resolution

  if (isTRUE(verbose) && n_steps > 0L) {
    message("EGA resolution search: explored ", n_steps + 1L,
            " resolution(s) from ", format(round(res_start, 3)),
            " downward; accepted resolution ",
            format(round(accepted_resolution, 3)),
            " (fewest items in communities below min_community_size = ",
            min_size, ").")
  }

  # Orphan any remaining sub-threshold communities (disconnected
  # components, or leftovers at the resolution floor).
  ega_native_na <- names(wc)[is.na(wc)]
  sizes <- table(wc[!is.na(wc)])
  small <- as.integer(names(sizes)[sizes < min_size])
  rule_na <- character(0)
  if (min_size > 1L && length(small) > 0) {
    rule_na <- names(wc)[!is.na(wc) & wc %in% small]
    wc[wc %in% small] <- NA_integer_
  }

  # Relabel surviving communities consecutively (1, 2, ...).
  if (any(!is.na(wc))) {
    wc[!is.na(wc)] <- as.integer(factor(wc[!is.na(wc)],
                                        levels = sort(unique(wc[!is.na(wc)]))))
  }
  n_dim <- length(unique(wc[!is.na(wc)]))

  if (all(is.na(wc))) {
    warning(
      "An empty item network was detected. This means that the input data likely",
      " do not contain a detectable common-factor or dimensional structure ",
      "(for example, random or noise-like item responses). EGA returned ",
      "0 dimensions and will not be used as a partition source unless ",
      "requested explicitly.",
      call. = FALSE
    )
    n_dim <- 0L
  } else if (anyNA(wc)) {
    parts <- character(0)
    if (length(ega_native_na) > 0) {
      parts <- c(parts, paste0(length(ega_native_na), " item(s) left ",
                               "unassigned by EGA itself (",
                               paste(ega_native_na, collapse = ", "), ")"))
    }
    if (length(rule_na) > 0) {
      parts <- c(parts, paste0(length(rule_na), " item(s) orphaned by the ",
                               "minimum-community-size rule (min_community_size = ",
                               min_size, "): ",
                               paste(rule_na, collapse = ", ")))
    }
    warning("EGA partition contains unassigned items: ",
            paste(parts, collapse = "; "),
            ". These items are excluded from the proposed measurement ",
            "syntax and SEM schematic.", call. = FALSE)
  }

  list(
    raw = ega,
    wc = wc,
    partition = wc,
    partition_table = .dimsem_partition_table(colnames(data), wc),
    n_dim = n_dim,
    network = network,
    resolution = accepted_resolution,
    resolution_steps = n_steps,
    community_sizes = table(wc[!is.na(wc)]),
    orphaned_small = rule_na
  )
}

.dimsem_run_pa <- function(data, pa_args = list(), seed = NULL) {
  .dimsem_require_namespace("psych")

  if (!is.null(seed)) {
    set.seed(seed)
  }

  defaults <- list(
    fm = "minres",
    fa = "fa",
    n.iter = 100,
    quant = 0.95,
    cor = "cor",
    use = "pairwise",
    plot = FALSE
  )

  args <- .dimsem_merge_args(defaults, pa_args)
  n_cores <- args$n_cores
  args$n_cores <- NULL

  # Convenience argument (intercepted, not passed to fa.parallel):
  # fa.parallel() deliberately overextracts factors from the observed and
  # resampled correlation matrices, which routinely produces (ultra-)Heywood
  # cases in those throwaway EFA solutions on discrete item data. The PA
  # retention decision compares eigenvalues and is unaffected, so these
  # warnings are muffled by default. Set pa_args$suppress_heywood = FALSE
  # to surface them; all other fa.parallel warnings always pass through.
  suppress_heywood <- args$suppress_heywood %||% TRUE
  args$suppress_heywood <- NULL

  old_mc <- getOption("mc.cores")
  on.exit(options(mc.cores = old_mc), add = TRUE)
  if (!is.null(n_cores)) {
    options(mc.cores = as.integer(n_cores))
  }

  run_pa <- function() do.call(psych::fa.parallel, c(list(x = data), args))
  pa <- if (isTRUE(suppress_heywood)) {
    withCallingHandlers(
      run_pa(),
      warning = function(w) {
        if (grepl("Heywood case", conditionMessage(w), fixed = TRUE)) {
          invokeRestart("muffleWarning")
        }
      }
    )
  } else {
    run_pa()
  }
  n_dim <- .dimsem_extract_pa_n_dim(pa, args)

  list(
    raw = pa,
    n_dim = .dimsem_clamp_dim(n_dim, ncol(data)),
    partition = NULL,
    partition_table = NULL,
    communities = NULL,
    network_matrix = NULL
  )
}

# Empirical Kaiser Criterion (Braeken & van Assen, 2017, Psychological
# Methods). Sample eigenvalues of the item correlation matrix are compared
# against serially updated reference eigenvalues derived from the
# Marchenko-Pastur upper bound, (1 + sqrt(p/n))^2, rescaled by the variance
# left unexplained by previously retained eigenvalues and floored at the
# classical Kaiser value of 1. Dimensions are retained while the sample
# eigenvalue exceeds its reference; counting stops at the first failure.
# Deterministic, closed form, base R only.
.dimsem_run_ekc <- function(data, ekc_args = list(), seed = NULL) {
  defaults <- list(use = "pairwise.complete.obs")
  args <- .dimsem_merge_args(defaults, ekc_args)

  n <- nrow(data)
  p <- ncol(data)
  R <- stats::cor(data, use = args$use)
  R[!is.finite(R)] <- 0
  diag(R) <- 1

  ev <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  ref <- numeric(p)
  for (j in seq_len(p)) {
    ref[j] <- max(
      ((1 + sqrt(p / n))^2) * (p - sum(ev[seq_len(j - 1)])) / (p - j + 1),
      1
    )
  }

  n_dim_raw <- sum(cumprod(ev > ref))

  list(
    raw = list(eigenvalues = ev, reference = ref),
    n_dim_raw = as.integer(n_dim_raw),
    # n_dim_raw = 0 (no retained dimensions, e.g. pure-noise data) is
    # clamped to a single dimension so a partition can still be built.
    n_dim = .dimsem_clamp_dim(max(n_dim_raw, 1L), p),
    partition = NULL,
    partition_table = NULL,
    communities = NULL,
    network_matrix = NULL
  )
}

# --- network estimation & graph partitions ---------------------------------------

.dimsem_get_network_matrix <- function(data,
                                       results,
                                       item_names,
                                       assignment_args = list()) {
  ega_network <- NULL
  if (!is.null(results$EGA$network)) {
    ega_network <- results$EGA$network
  }

  if (!is.null(ega_network)) {
    return(.dimsem_named_square_matrix(ega_network, item_names))
  }

  .dimsem_estimate_glasso_network(
    data = data,
    item_names = item_names,
    assignment_args = assignment_args
  )
}

.dimsem_estimate_glasso_network <- function(data,
                                            item_names,
                                            assignment_args = list()) {
  .dimsem_require_namespace("qgraph")

  glasso_args <- assignment_args$glasso_args %||% list()
  cor_method <- assignment_args$cor %||% "auto"
  use <- assignment_args$use %||% "pairwise.complete.obs"

  if (is.matrix(cor_method) || is.data.frame(cor_method)) {
    cor_mat <- as.matrix(cor_method)
  } else if (identical(tolower(cor_method), "auto")) {
    cor_mat <- qgraph::cor_auto(data)
  } else if (identical(tolower(cor_method), "pearson")) {
    cor_mat <- stats::cor(data, use = use)
  } else {
    stop("`assignment_args$cor` must be \"auto\", \"pearson\", or a matrix.",
         call. = FALSE)
  }

  cor_mat <- .dimsem_named_square_matrix(cor_mat, item_names)
  defaults <- list(gamma = 0.5)
  glasso_args <- .dimsem_merge_args(defaults, glasso_args)
  ebic_args <- .dimsem_merge_args(
    list(S = cor_mat, n = nrow(data)),
    glasso_args
  )
  network <- do.call(qgraph::EBICglasso, ebic_args)
  .dimsem_named_square_matrix(network, item_names)
}

# Attach a k-community item partition to every method that only estimates a
# dimension *count* (PA, EKC). EGA supplies its own communities.
.dimsem_add_graph_partitions <- function(results,
                                         network_matrix,
                                         item_names,
                                         assignment_args = list(),
                                         seed = NULL,
                                         verbose = TRUE) {
  algorithm <- assignment_args$algorithm %||% "walktrap"
  needs_partition <- list(
    PA = results$PA$n_dim,
    EKC = results$EKC$n_dim
  )

  for (method in names(needs_partition)) {
    n_dim <- needs_partition[[method]]
    if (is.null(results[[method]]) || is.null(n_dim)) {
      next
    }
    if (isTRUE(verbose)) {
      message("Assigning ", method, " dimensions with weighted graph ",
              "partitioning (", algorithm, ")...")
    }
    part <- .dimsem_k_partition(
      network_matrix = network_matrix,
      n_dim = n_dim,
      item_names = item_names,
      assignment_args = assignment_args,
      seed = seed
    )
    results[[method]]$partition <- part$partition
    results[[method]]$partition_table <- part$partition_table
    results[[method]]$communities <- part$communities
    results[[method]]$network_matrix <- network_matrix
    results[[method]]$n_dim_partition <- length(unique(part$partition))
    results[[method]]$partition_algorithm <- part$algorithm
  }

  results
}

# Partition the item network into exactly `n_dim` communities.
#
# Algorithms (assignment_args$algorithm):
#   "walktrap" (default) : igraph::cluster_walktrap() with edge weights,
#                          cut at exactly n_dim communities via
#                          igraph::cut_at(). Deterministic and
#                          weight-respecting.
#   "fast_greedy"        : igraph::cluster_fast_greedy() with edge weights,
#                          cut at n_dim via igraph::cut_at(). Deterministic.
#   "fluid"              : igraph::cluster_fluid_communities(). NOTE: this
#                          algorithm ignores edge weights entirely; it is
#                          retained for backward compatibility but is not
#                          recommended for weighted regularized networks.
.dimsem_k_partition <- function(network_matrix,
                                n_dim,
                                item_names,
                                assignment_args = list(),
                                seed = NULL) {
  network_matrix <- .dimsem_named_square_matrix(network_matrix, item_names)
  p <- ncol(network_matrix)
  n_dim <- .dimsem_clamp_dim(n_dim, p)

  algorithm <- assignment_args$algorithm %||% "walktrap"
  algorithm <- match.arg(algorithm,
                         choices = c("walktrap", "fast_greedy", "fluid"))

  if (n_dim <= 1) {
    partition <- stats::setNames(rep(1L, p), item_names)
    return(list(
      partition = partition,
      partition_table = .dimsem_partition_table(item_names, partition),
      communities = NULL,
      graph = NULL,
      algorithm = algorithm
    ))
  }

  if (n_dim >= p) {
    partition <- stats::setNames(seq_len(p), item_names)
    return(list(
      partition = partition,
      partition_table = .dimsem_partition_table(item_names, partition),
      communities = NULL,
      graph = NULL,
      algorithm = algorithm
    ))
  }

  g <- .dimsem_partition_graph(network_matrix, assignment_args = assignment_args)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  detect <- function() {
    if (identical(algorithm, "fluid")) {
      # Backward-compatible option. cluster_fluid_communities() ignores edge
      # weights, so strong within-dimension and weak cross-dimension edges
      # count equally; prefer "walktrap" for weighted regularized networks.
      communities <- igraph::cluster_fluid_communities(g, n_dim)
      return(list(communities = communities,
                  membership = igraph::membership(communities)))
    }

    communities <- switch(
      algorithm,
      walktrap = igraph::cluster_walktrap(g, weights = igraph::E(g)$weight),
      fast_greedy = igraph::cluster_fast_greedy(g,
                                                weights = igraph::E(g)$weight)
    )
    # Both algorithms are agglomerative; cut their merge dendrogram at
    # exactly n_dim communities.
    list(communities = communities,
         membership = igraph::cut_at(communities, no = n_dim))
  }

  detected <- tryCatch(
    detect(),
    error = function(e) {
      stop("Community detection (", algorithm, ") failed: ",
           conditionMessage(e), call. = FALSE)
    }
  )

  partition <- detected$membership
  partition <- as.integer(factor(partition, levels = sort(unique(partition))))
  names(partition) <- igraph::V(g)$name
  partition <- partition[item_names]

  list(
    partition = partition,
    partition_table = .dimsem_partition_table(item_names, partition),
    communities = detected$communities,
    graph = g,
    algorithm = algorithm
  )
}

.dimsem_partition_graph <- function(network_matrix, assignment_args = list()) {
  transform <- assignment_args$weight_transform %||% "absolute"
  transform <- match.arg(transform, choices = c("absolute", "positive"))
  connect <- assignment_args$connect_disconnected %||% TRUE

  weights <- network_matrix
  weights[!is.finite(weights)] <- 0
  diag(weights) <- 0
  if (identical(transform, "absolute")) {
    weights <- abs(weights)
  } else {
    weights <- pmax(weights, 0)
  }

  weights[weights < .Machine$double.eps] <- 0
  g <- igraph::graph_from_adjacency_matrix(
    weights,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )
  g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)

  if (igraph::ecount(g) == 0) {
    g <- igraph::make_full_graph(igraph::vcount(g), directed = FALSE)
    igraph::V(g)$name <- colnames(network_matrix)
    igraph::E(g)$weight <- .Machine$double.eps
    return(g)
  }

  if (isTRUE(connect) && !igraph::is_connected(g)) {
    g <- .dimsem_connect_components(g)
  }

  g
}

.dimsem_connect_components <- function(g) {
  comps <- igraph::components(g)
  if (comps$no <= 1) {
    return(g)
  }

  reps <- vapply(
    split(seq_along(comps$membership), comps$membership),
    function(idx) idx[1],
    integer(1)
  )
  edges <- as.vector(rbind(reps[-length(reps)], reps[-1]))
  existing_weights <- igraph::E(g)$weight
  positive_weights <- existing_weights[is.finite(existing_weights) &
                                         existing_weights > 0]
  bridge_weight <- if (length(positive_weights) > 0) {
    min(positive_weights) * 1e-6
  } else {
    .Machine$double.eps
  }

  old_ecount <- igraph::ecount(g)
  g <- igraph::add_edges(g, edges)
  new_edges <- seq.int(old_ecount + 1L, igraph::ecount(g))
  igraph::E(g)$weight[new_edges] <- bridge_weight
  g
}

# --- partition selection ---------------------------------------------------------

.dimsem_select_partition_source <- function(requested, results) {
  available <- names(results)

  has_partition <- vapply(
    results,
    function(res) !is.null(res$partition),
    logical(1)
  )

  usable <- vapply(
    results,
    function(res) !is.null(res$partition) && any(!is.na(res$partition)),
    logical(1)
  )

  empty_ega <- (
    "EGA" %in% available &&
      isTRUE(has_partition[["EGA"]]) &&
      !is.null(results$EGA$n_dim) &&
      identical(as.integer(results$EGA$n_dim), 0L) &&
      all(is.na(results$EGA$partition))
  )

  if (!identical(requested, "auto")) {
    if (!requested %in% available) {
      stop("Requested `partition_source = \"", requested,
           "\"`, but that method was not run. Available results: ",
           paste(available, collapse = ", "), call. = FALSE)
    }

    if (usable[[requested]]) {
      return(requested)
    }

    if (identical(requested, "EGA") && isTRUE(empty_ega)) {
      return("EGA")
    }

    stop("Requested `partition_source = \"", requested, "\"`, but that ",
         "method did not produce a usable item partition. Usable sources: ",
         if (any(usable)) paste(available[usable], collapse = ", ")
         else "none",
         ".", call. = FALSE)
  }

  for (candidate in c("EGA", "PA", "EKC")) {
    if (candidate %in% available && usable[[candidate]]) {
      return(candidate)
    }
  }

  if (isTRUE(empty_ega)) {
    return("EGA")
  }

  stop("No method produced a usable item partition.", call. = FALSE)
}

.dimsem_build_selected_partition <- function(item_names,
                                             results,
                                             selected_source) {
  result <- results[[selected_source]]
  partition <- result$partition
  if (is.null(partition)) {
    stop("Selected method does not contain an item partition.", call. = FALSE)
  }

  partition <- partition[item_names]
  partition_table <- .dimsem_partition_table(item_names, partition)
  n_dim <- length(unique(partition[!is.na(partition)]))

  list(
    source = selected_source,
    n_dim = n_dim,
    partition = partition,
    partition_table = partition_table,
    communities = result$communities %||% NULL,
    model_syntax = .dimsem_model_syntax(partition_table)
  )
}

.dimsem_partition_table <- function(item_names, partition) {
  data.frame(
    item = item_names,
    dimension = as.integer(partition[item_names]),
    stringsAsFactors = FALSE
  )
}

.dimsem_model_syntax <- function(partition_table) {
  tab <- partition_table[!is.na(partition_table$dimension), , drop = FALSE]
  if (nrow(tab) == 0L) {
    return("")
  }
  dims <- sort(unique(tab$dimension))
  lines <- vapply(
    dims,
    function(d) {
      items <- tab$item[tab$dimension == d]
      paste0("F", d, " =~ ", paste(items, collapse = " + "))
    },
    character(1)
  )

  paste(lines, collapse = "\n")
}

# --- confirmatory factor analysis -------------------------------------------------

# Fit a lavaan CFA on one of the proposed item partitions and extract the
# inter-factor correlation matrix.
#
#   cfa_source : "selected" (default) uses the partition already chosen via
#                `partition_source` (EGA under auto-selection whenever EGA ran
#                and produced communities); "EGA", "PA", or "EKC" force the
#                CFA onto that method's partition instead.
#   cfa_args   : named list merged over `list(std.lv = TRUE)` and passed to
#                lavaan::cfa(). std.lv = TRUE standardizes the latent
#                variables so latent covariances are correlations and every
#                loading is freely estimated. For ordinal treatment of
#                Likert items, pass list(ordered = TRUE) (lavaan then uses
#                polychoric correlations with a WLSMV-type estimator).
.dimsem_fit_cfa <- function(data,
                            results,
                            selected,
                            cfa_source = "selected",
                            cfa_args = list()) {
  .dimsem_require_namespace("lavaan")

  if (identical(cfa_source, "selected")) {
    source_used <- selected$source
    partition_table <- selected$partition_table
    syntax <- selected$model_syntax
  } else {
    if (!cfa_source %in% names(results)) {
      stop("Requested `cfa_source = \"", cfa_source, "\"`, but that method ",
           "was not run. Available results: ",
           paste(names(results), collapse = ", "), call. = FALSE)
    }
    res <- results[[cfa_source]]
    if (is.null(res$partition) || all(is.na(res$partition))) {
      stop("Requested `cfa_source = \"", cfa_source, "\"`, but that method ",
           "did not produce a usable item partition.", call. = FALSE)
    }
    source_used <- cfa_source
    partition_table <- res$partition_table
    syntax <- .dimsem_model_syntax(partition_table)
  }

  cfa_out <- list(
    source = source_used,
    model_syntax = syntax,
    fit = NULL,
    converged = FALSE,
    factor_cor = NULL,
    fit_measures = NULL,
    error = NULL
  )

  if (!nzchar(syntax)) {
    warning("No CFA was fitted: the ", source_used, " partition assigns no ",
            "items to any dimension (empty measurement model).",
            call. = FALSE)
    cfa_out$error <- "empty measurement model"
    return(cfa_out)
  }

  # Identification heads-up (lavaan may still fit with equality constraints
  # borrowed across factors, but users should look at such factors closely).
  assigned <- partition_table[!is.na(partition_table$dimension), , drop = FALSE]
  dim_sizes <- table(assigned$dimension)
  n_factors <- length(dim_sizes)
  if (n_factors == 1L && dim_sizes[[1]] < 3L) {
    warning("The single proposed factor has fewer than three indicators; ",
            "the CFA is not identified.", call. = FALSE)
  } else if (any(dim_sizes < 2L)) {
    warning("Dimension(s) ", paste(names(dim_sizes)[dim_sizes < 2L],
                                   collapse = ", "), " have only one indicator; interpret the ",
            "CFA solution with caution.", call. = FALSE)
  }

  # ordered = TRUE by default: dimsem's target data are Likert-type ordered
  # (incl. binary) items, for which lavaan then uses polychoric/tetrachoric
  # correlations with a WLSMV-type (DWLS + robust corrections) estimator.
  # Override with cfa_args = list(ordered = FALSE) to treat items as
  # continuous.
  defaults <- list(std.lv = TRUE, ordered = TRUE, missing = "pairwise")
  args <- .dimsem_merge_args(defaults, cfa_args)

  # lavaan::cfa() re-evaluates parts of its captured call (non-standard
  # evaluation); splicing the data.frame *value* in via do.call() makes that
  # re-evaluation resolve `data` to the utils::data() closure and fail with
  # "cannot coerce type 'closure'". Building the call with a quoted symbol
  # that is looked up in this function's frame avoids the problem.
  lav_call <- as.call(c(
    list(quote(lavaan::cfa)),
    list(model = syntax, data = quote(data)),
    args
  ))
  fit <- tryCatch(eval(lav_call), error = function(e) e)

  if (inherits(fit, "error")) {
    warning("lavaan CFA estimation failed: ", conditionMessage(fit),
            call. = FALSE)
    cfa_out$error <- conditionMessage(fit)
    return(cfa_out)
  }

  cfa_out$fit <- fit
  cfa_out$converged <- isTRUE(tryCatch(
    lavaan::lavInspect(fit, "converged"),
    error = function(e) FALSE
  ))

  if (!cfa_out$converged) {
    warning("The lavaan CFA did not converge; factor correlations and fit ",
            "measures are not reported. The (non-converged) fit object is ",
            "stored in `$cfa$fit` for inspection.", call. = FALSE)
    return(cfa_out)
  }

  # With std.lv = TRUE, "cov.lv" already contains correlations; "cor.lv"
  # is robust to user-supplied overrides of std.lv as well.
  cfa_out$factor_cor <- tryCatch(
    lavaan::lavInspect(fit, "cor.lv"),
    error = function(e) NULL
  )
  cfa_out$fit_measures <- tryCatch({
    fm <- lavaan::fitMeasures(fit)
    plain <- c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr")
    # Under ordered/WLSMV estimation, the mean-and-variance-adjusted
    # (scaled/robust) statistics are the ones to report; fall back to the
    # plain statistic when no scaled variant exists (e.g., srmr, or ML fits).
    scaled <- paste0(plain, ".scaled")
    picked <- ifelse(scaled %in% names(fm) & !is.na(fm[scaled]), scaled, plain)
    out <- fm[picked]
    names(out) <- plain
    attr(out, "scaled") <- any(picked != plain)
    out
  },
  error = function(e) NULL
  )

  cfa_out
}

# --- graphs & plots ---------------------------------------------------------------

# Choose the hyper result used for plotting.
#
# Plotting priority rule:
#   1. permutation / higher-order Louvain whenever available
#   2. TGA only when permutation was not run
#
# This means that if TGA and permutation disagree, the plot follows the
# permutation/higher-order Louvain solution. This helper is plotting-only;
# it does not alter the hyper proposal object or recommendation logic.
.dimsem_hyper_plot_result <- function(hyper) {
  if (is.null(hyper)) {
    return(NULL)
  }

  if (!is.null(hyper$permutation)) {
    res <- hyper$permutation
    source <- "permutation"
  } else if (!is.null(hyper$tga)) {
    res <- hyper$tga
    source <- "tga"
  } else {
    return(NULL)
  }

  blocks <- res$blocks %||% list()

  if (length(blocks) > 0L) {
    empty <- vapply(blocks, length, integer(1)) == 0L
    blocks <- blocks[!empty]
  }

  if (length(blocks) > 0L &&
      (is.null(names(blocks)) || anyNA(names(blocks)) ||
       any(names(blocks) == ""))) {
    names(blocks) <- paste0("G", seq_along(blocks))
  }

  list(
    source = source,
    verdict = res$verdict %||% NA_character_,
    blocks = blocks,
    result = res
  )
}

.dimsem_factor_hyper_membership <- function(latent_names, hyper) {
  out <- stats::setNames(rep(NA_character_, length(latent_names)), latent_names)
  hp <- .dimsem_hyper_plot_result(hyper)

  if (is.null(hp) || length(hp$blocks) == 0L) {
    return(out)
  }

  for (g in names(hp$blocks)) {
    members <- intersect(hp$blocks[[g]], latent_names)
    out[members] <- g
  }

  out
}

.dimsem_item_hyper_membership <- function(partition, hyper) {
  out <- stats::setNames(rep(NA_character_, length(partition)), names(partition))
  hp <- .dimsem_hyper_plot_result(hyper)

  if (is.null(hp) || length(hp$blocks) == 0L) {
    return(out)
  }

  dim_label <- ifelse(
    is.na(partition),
    NA_character_,
    paste0("F", as.integer(partition))
  )

  for (g in names(hp$blocks)) {
    out[!is.na(dim_label) & dim_label %in% hp$blocks[[g]]] <- g
  }

  out
}

.dimsem_network_graph <- function(network_matrix, partition, hyper = NULL) {
  item_names <- colnames(network_matrix)
  network_matrix <- .dimsem_named_square_matrix(network_matrix, item_names)
  partition <- partition[item_names]

  signed <- network_matrix
  signed[!is.finite(signed)] <- 0
  diag(signed) <- 0

  g <- igraph::graph_from_adjacency_matrix(
    abs(signed),
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )
  g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)

  v_names <- igraph::V(g)$name
  v_dim <- as.integer(partition[v_names])

  igraph::V(g)$dimension <- v_dim
  igraph::V(g)$dimension_label <- ifelse(
    is.na(v_dim),
    NA_character_,
    paste0("F", v_dim)
  )

  item_hyper <- .dimsem_item_hyper_membership(partition, hyper)
  igraph::V(g)$hyper <- item_hyper[v_names]

  hp <- .dimsem_hyper_plot_result(hyper)
  igraph::graph_attr(g, "hyper_source") <- if (is.null(hp)) NA_character_ else hp$source
  igraph::graph_attr(g, "has_hyper") <- !is.null(hp) && length(hp$blocks) > 0L

  if (igraph::ecount(g) > 0) {
    edge_ends <- igraph::ends(g, igraph::E(g), names = TRUE)
    signed_weights <- signed[cbind(edge_ends[, 1], edge_ends[, 2])]
    igraph::E(g)$signed_weight <- signed_weights
    igraph::E(g)$abs_weight <- abs(signed_weights)
    igraph::E(g)$sign <- ifelse(signed_weights >= 0, "positive", "negative")
  }

  g
}

.dimsem_sem_graph <- function(partition_table, hyper = NULL) {
  assigned <- partition_table[!is.na(partition_table$dimension), , drop = FALSE]

  item_vertices <- data.frame(
    name = partition_table$item,
    type = "item",
    hyper = NA_character_,
    stringsAsFactors = FALSE
  )

  if (nrow(assigned) == 0L) {
    g <- igraph::graph_from_data_frame(
      d = data.frame(
        from = character(0),
        to = character(0),
        edge_type = character(0),
        stringsAsFactors = FALSE
      ),
      directed = TRUE,
      vertices = item_vertices
    )
    igraph::graph_attr(g, "has_hyper") <- FALSE
    igraph::graph_attr(g, "hyper_source") <- NA_character_
    return(g)
  }

  latent_names <- paste0("F", sort(unique(assigned$dimension)))
  latent_hyper <- .dimsem_factor_hyper_membership(latent_names, hyper)

  latent_vertices <- data.frame(
    name = latent_names,
    type = "latent",
    hyper = unname(latent_hyper[latent_names]),
    stringsAsFactors = FALSE
  )

  hp <- .dimsem_hyper_plot_result(hyper)

  hyper_vertices <- data.frame(
    name = character(0),
    type = character(0),
    hyper = character(0),
    stringsAsFactors = FALSE
  )

  hyper_edges <- data.frame(
    from = character(0),
    to = character(0),
    edge_type = character(0),
    stringsAsFactors = FALSE
  )

  if (!is.null(hp) && length(hp$blocks) > 0L) {
    hyper_names <- names(hp$blocks)

    hyper_vertices <- data.frame(
      name = hyper_names,
      type = "hyper",
      hyper = hyper_names,
      stringsAsFactors = FALSE
    )

    hyper_edges <- do.call(
      rbind,
      lapply(hyper_names, function(g) {
        children <- intersect(hp$blocks[[g]], latent_names)
        if (length(children) == 0L) {
          return(NULL)
        }
        data.frame(
          from = g,
          to = children,
          edge_type = "higher_order",
          stringsAsFactors = FALSE
        )
      })
    )

    if (is.null(hyper_edges)) {
      hyper_edges <- data.frame(
        from = character(0),
        to = character(0),
        edge_type = character(0),
        stringsAsFactors = FALSE
      )
    }
  }

  loading_edges <- data.frame(
    from = paste0("F", assigned$dimension),
    to = assigned$item,
    edge_type = "loading",
    stringsAsFactors = FALSE
  )

  vertices <- rbind(hyper_vertices, latent_vertices, item_vertices)
  edges <- rbind(hyper_edges, loading_edges)

  g <- igraph::graph_from_data_frame(edges, directed = TRUE, vertices = vertices)

  igraph::graph_attr(g, "has_hyper") <- !is.null(hp) && length(hp$blocks) > 0L
  igraph::graph_attr(g, "hyper_source") <- if (is.null(hp)) NA_character_ else hp$source

  g
}

.dimsem_plot_network <- function(graph, plot_args = list()) {
  if (!requireNamespace("ggraph", quietly = TRUE) ||
      !requireNamespace("ggplot2", quietly = TRUE)) {
    warning("Packages `ggraph` and `ggplot2` are required for plotting.",
            call. = FALSE)
    return(NULL)
  }

  layout <- plot_args$network_layout %||% "fr"
  node_size <- plot_args$node_size %||% 4
  label_size <- plot_args$label_size %||% 3

  edge_colours <- plot_args$edge_colours %||% c(
    positive = "grey35",
    negative = "grey75"
  )

  show_hyper_hulls <- plot_args$hyper_hulls %||% TRUE
  hyper_alpha <- plot_args$hyper_alpha %||% 0.14
  hyper_expand <- plot_args$hyper_expand %||% 4
  hyper_radius <- plot_args$hyper_radius %||% 4
  hyper_concavity <- plot_args$hyper_concavity %||% 5
  hyper_min_points <- plot_args$hyper_min_points %||% 3L

  lay <- ggraph::create_layout(graph, layout = layout)
  p <- ggraph::ggraph(lay)

  if (isTRUE(show_hyper_hulls) &&
      "hyper" %in% names(as.data.frame(lay))) {
    hull_data <- as.data.frame(lay)
    hull_data <- hull_data[
      !is.na(hull_data$hyper) & hull_data$hyper != "",
      ,
      drop = FALSE
    ]

    if (nrow(hull_data) > 0L) {
      counts <- table(hull_data$hyper)
      keep <- names(counts)[counts >= hyper_min_points]
      hull_data <- hull_data[hull_data$hyper %in% keep, , drop = FALSE]
    }

    if (nrow(hull_data) > 0L) {
      if (requireNamespace("ggforce", quietly = TRUE)) {
        p <- p +
          ggforce::geom_mark_hull(
            data = hull_data,
            ggplot2::aes(
              x = .data$x,
              y = .data$y,
              group = .data$hyper,
              fill = .data$hyper,
              label = .data$hyper
            ),
            concavity = hyper_concavity,
            expand = grid::unit(hyper_expand, "mm"),
            radius = grid::unit(hyper_radius, "mm"),
            alpha = hyper_alpha,
            show.legend = FALSE,
            inherit.aes = FALSE
          )
      } else {
        warning(
          "Package `ggforce` is required for hypergraph-style background ",
          "hulls. Install `ggforce` or set `plot_args = list(hyper_hulls = FALSE)`.",
          call. = FALSE
        )
      }
    }
  }

  if (igraph::ecount(graph) > 0) {
    p <- p +
      ggraph::geom_edge_link(
        ggplot2::aes(
          edge_alpha = .data$abs_weight,
          edge_width = .data$abs_weight,
          edge_colour = factor(.data$sign)
        ),
        show.legend = FALSE
      ) +
      ggraph::scale_edge_colour_manual(
        values = edge_colours,
        drop = FALSE
      )
  }

  p +
    ggraph::geom_node_point(
      ggplot2::aes(color = factor(.data$dimension)),
      size = node_size
    ) +
    ggraph::geom_node_text(
      ggplot2::aes(label = .data$name),
      repel = TRUE,
      size = label_size
    ) +
    ggplot2::labs(
      color = "Dimension",
      fill = "Hyper-dimension"
    ) +
    ggplot2::theme_void()
}

.dimsem_plot_sem <- function(graph, plot_args = list()) {
  if (!requireNamespace("ggraph", quietly = TRUE) ||
      !requireNamespace("ggplot2", quietly = TRUE)) {
    warning("Packages `ggraph` and `ggplot2` are required for plotting.",
            call. = FALSE)
    return(NULL)
  }

  layout <- plot_args$sem_layout %||% "sugiyama"
  node_size <- plot_args$node_size %||% 5
  label_size <- plot_args$label_size %||% 3

  ggraph::ggraph(graph, layout = layout) +
    ggraph::geom_edge_link(
      arrow = grid::arrow(length = grid::unit(3, "mm"), type = "closed"),
      end_cap = ggraph::circle(3, "mm"),
      show.legend = FALSE
    ) +
    ggraph::geom_node_point(
      ggplot2::aes(shape = .data$type),
      size = node_size
    ) +
    ggraph::geom_node_text(
      ggplot2::aes(label = .data$name),
      repel = TRUE,
      size = label_size
    ) +
    ggplot2::labs(shape = "Node type") +
    ggplot2::theme_void()
}

# --- result extractors -------------------------------------------------------------

.dimsem_extract_ega_wc <- function(ega) {
  candidates <- list(
    ega$wc,
    ega$community,
    ega$communities,
    ega$dimension,
    ega$dim.variables
  )
  candidates <- candidates[vapply(candidates, function(x) !is.null(x), logical(1))]
  if (length(candidates) == 0) {
    stop("Could not find an EGA community vector in the returned object.",
         call. = FALSE)
  }
  wc <- candidates[[1]]
  as.integer(wc)
}

.dimsem_extract_ega_n_dim <- function(ega, wc) {
  candidates <- list(ega$n.dim, ega$n_dim, ega$ndim, ega$dimensions)
  candidates <- candidates[vapply(candidates, function(x) !is.null(x), logical(1))]
  if (length(candidates) > 0 && length(candidates[[1]]) == 1) {
    return(as.integer(candidates[[1]]))
  }

  length(unique(wc))
}

.dimsem_extract_ega_network <- function(ega, item_names) {
  preferred <- list(
    ega$network,
    ega$network.matrix,
    ega$glasso,
    ega$EGA$network
  )
  for (candidate in preferred) {
    candidate <- .dimsem_find_square_matrix(candidate, length(item_names))
    if (!is.null(candidate)) {
      return(.dimsem_named_square_matrix(candidate, item_names))
    }
  }

  .dimsem_find_square_matrix(ega, length(item_names))
}

.dimsem_find_square_matrix <- function(x, p, depth = 0) {
  if (is.null(x) || depth > 3) {
    return(NULL)
  }
  if (is.matrix(x) || is.data.frame(x)) {
    x <- as.matrix(x)
    if (nrow(x) == p && ncol(x) == p) {
      return(x)
    }
    return(NULL)
  }
  if (!is.list(x)) {
    return(NULL)
  }

  for (element in x) {
    candidate <- .dimsem_find_square_matrix(element, p, depth = depth + 1L)
    if (!is.null(candidate)) {
      return(candidate)
    }
  }

  NULL
}

.dimsem_extract_pa_n_dim <- function(pa, args) {
  fa_mode <- args$fa %||% "fa"
  candidates <- if (identical(fa_mode, "pc")) {
    list(pa$ncomp, pa$n.comp, pa$components, pa$nfact)
  } else {
    list(pa$nfact, pa$nfactors, pa$factors, pa$ncomp)
  }

  candidates <- candidates[vapply(
    candidates,
    function(x) !is.null(x) && length(x) >= 1 && is.finite(x[1]),
    logical(1)
  )]
  if (length(candidates) == 0) {
    stop("Could not extract a PA dimension count from `psych::fa.parallel()`.",
         call. = FALSE)
  }

  as.integer(round(candidates[[1]][1]))
}

# --- small utilities -----------------------------------------------------------------

.dimsem_named_square_matrix <- function(x, item_names) {
  x <- as.matrix(x)
  if (nrow(x) != length(item_names) || ncol(x) != length(item_names)) {
    stop("Network matrix must be square with one row/column per item.",
         call. = FALSE)
  }
  storage.mode(x) <- "numeric"
  rownames(x) <- item_names
  colnames(x) <- item_names
  x
}

.dimsem_clamp_dim <- function(n_dim, p) {
  if (length(n_dim) != 1 || !is.finite(n_dim)) {
    stop("Dimension count must be a finite scalar.", call. = FALSE)
  }

  max(1L, min(as.integer(ceiling(n_dim)), as.integer(p)))
}

.dimsem_merge_args <- function(defaults, user) {
  if (is.null(user)) user <- list()
  if (!is.list(user)) {
    stop("Argument lists must be named lists.", call. = FALSE)
  }
  utils::modifyList(defaults, user, keep.null = TRUE)
}

.dimsem_require_namespace <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package `", pkg, "` is required for this operation. ",
         "Please install it first.", call. = FALSE)
  }
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

# Resolve a requested core count into a usable one. Defaults to
# parallel::detectCores() - 1 (leaving one core for background tasks),
# clamps to [1, available], and never requests more workers than there are
# tasks to run.
.dimsem_resolve_ncores <- function(ncores = NULL, n_tasks = Inf) {
  if (!requireNamespace("parallel", quietly = TRUE)) {
    return(1L)
  }
  avail <- tryCatch(parallel::detectCores(), error = function(e) NA_integer_)
  if (!is.finite(avail) || avail < 1L) {
    avail <- 1L
  }
  if (is.null(ncores)) {
    ncores <- avail - 1L
  }
  ncores <- suppressWarnings(as.integer(ncores)[1])
  if (!is.finite(ncores) || ncores < 1L) {
    ncores <- 1L
  }
  min(ncores, avail, max(1L, n_tasks))
}

# Hyper-dimensional (level-2) structure proposal for dimsem ---------------------
#
# Two network-psychometric strategies for deciding whether the level-1 item
# partition warrants one or more "hyper-factors" (second-order dimensions):
#
#   (1) "tga": a two-level adaptation of Taxonomic Graph Analysis (Samo,
#       Garrido, Abad, Golino, McAbee, & Christensen, 2025, European Journal
#       of Personality). Network loadings (Christensen, Golino, Abad, &
#       Garrido, 2025) are computed from the level-1 network and partition,
#       network scores (Golino et al., 2022) are correlated, a level-2
#       network is estimated on the scores, and lower-order Louvain community
#       detection (Blondel et al., 2008; Jimenez et al., 2023) with
#       "most common" consensus clustering (Lancichinetti & Fortunato, 2012;
#       Golino & Christensen, 2025) identifies level-2 communities. Whether a
#       single top-level dimension is statistically supported is evaluated
#       with the unidim index (Revelle & Condon, 2025).
#       NOTE: this borrows TGA's hierarchical core only. Full TGA includes
#       item-hygiene steps (UVA for local dependence, riEGA for wording
#       effects) and bootstrap robustness checks that are out of scope here;
#       users can pre-clean items with EGAnet before calling dimsem_propose().
#
#   (2) "permutation": a statistically validated higher-order Louvain. The
#       null hypothesis is that level-1 communities are MUTUALLY INDEPENDENT
#       conditional on their internal structure. This null is simulated
#       exactly and nonparametrically by permuting data rows independently
#       within each community's column block, which preserves every
#       community's internal joint distribution while destroying all
#       between-community dependence. The test statistic for each community
#       pair is the between-community edge mass of the re-estimated
#       regularized network, S_cd = sum_{i in c, j in d} |w_ij| -- i.e., the
#       between-community component of the weighted-modularity numerator.
#       Pairs whose observed S_cd exceeds the permutation null (BH-FDR
#       corrected) form the weighted edges of a meta-graph over communities;
#       Louvain consensus on this meta-graph yields the hyper-blocks, and
#       communities with no validated edges remain stand-alone ("orphans").
#       A global permutation p-value on the total between-community edge
#       mass acts as an omnibus gate.

# --- network loadings & scores --------------------------------------------------

# Signed, standardized network loadings (Christensen & Golino, 2021;
# Christensen et al., 2025). Unstandardized loading of node i on community c
# is the sum of its (absolute) edge weights to the members of c; the sign is
# taken from the signed sum. Standardization divides by the square root of
# the total within-community loading mass of c's own members, which puts
# loadings on a scale comparable to factor loadings.
.dimsem_network_loadings <- function(network_matrix, partition) {
  part <- partition[!is.na(partition)]
  comms <- sort(unique(part))
  items <- names(part)
  W <- network_matrix[items, items, drop = FALSE]
  diag(W) <- 0

  L_abs <- sapply(comms, function(cc) {
    members <- items[part == cc]
    rowSums(abs(W[, members, drop = FALSE]))
  })
  L_sgn <- sapply(comms, function(cc) {
    members <- items[part == cc]
    rowSums(W[, members, drop = FALSE])
  })
  L_abs <- matrix(L_abs, nrow = length(items),
                  dimnames = list(items, paste0("F", comms)))
  L_sgn <- matrix(L_sgn, nrow = length(items),
                  dimnames = list(items, paste0("F", comms)))

  L_std <- L_abs
  for (k in seq_along(comms)) {
    members <- items[part == comms[k]]
    denom <- sqrt(sum(L_abs[members, k]))
    L_std[, k] <- if (denom > 0) L_abs[, k] / denom else 0
  }

  sgn <- sign(L_sgn)
  sgn[sgn == 0] <- 1
  L_std * sgn
}

# Network scores (Golino et al., 2022): standardized data weighted by each
# item's loading on its ASSIGNED community only, then re-standardized.
.dimsem_network_scores <- function(data, network_matrix, partition) {
  part <- partition[!is.na(partition)]
  comms <- sort(unique(part))
  items <- names(part)

  L <- .dimsem_network_loadings(network_matrix, partition)
  Z <- scale(as.matrix(data[, items, drop = FALSE]))
  # scale() maps a zero-variance column to NaN; those carry no information
  # and are set to the (standardized) mean of 0. Genuine NAs are NOT
  # touched here -- they are handled per-respondent below, so that
  # incomplete rows are neither deleted nor silently mean-imputed.
  zero_var <- vapply(seq_len(ncol(Z)), function(j) {
    all(is.na(Z[, j])) || (stats::sd(data[[items[j]]], na.rm = TRUE) %in% c(0, NA))
  }, logical(1))
  if (any(zero_var)) {
    Z[, zero_var] <- 0
  }

  scores <- sapply(seq_along(comms), function(k) {
    members <- items[part == comms[k]]
    Zm <- Z[, members, drop = FALSE]
    w <- L[members, k]
    obs <- !is.na(Zm)
    # Weighted sum over observed items, rescaled by the share of loading
    # mass observed for that respondent (|w| in the denominator keeps the
    # rescaling sign-safe for negatively loading items).
    num <- as.numeric(replace(Zm, !obs, 0) %*% w)
    den <- as.numeric(obs %*% abs(w))
    total <- sum(abs(w))
    out <- ifelse(den > 0, num * (total / pmax(den, .Machine$double.eps)), NA_real_)
    out
  })
  scores <- matrix(scores, ncol = length(comms),
                   dimnames = list(NULL, paste0("F", comms)))
  # Standardize using observed values only; constant/all-NA columns stay 0.
  scores <- scale(scores)
  const <- vapply(seq_len(ncol(scores)), function(j) all(!is.finite(scores[, j])),
                  logical(1))
  if (any(const)) {
    scores[, const] <- 0
  }
  scores
}

# --- unidimensionality index ------------------------------------------------------

# Internal implementation of the unidim index (Revelle & Condon, 2025):
# u = tau * pc, the product of a tau-equivalence fit index (single factor,
# equal loadings) and a congeneric fit index (single factor, free loadings),
# both computed on the off-diagonal correlation structure. Values near 1
# support unidimensionality; Revelle & Condon's simulations place strong
# support around >= .90 and clear multidimensionality around <= .50.
# psych::unidim() is the reference implementation; this internal version
# avoids cross-version fragility and requires only psych::fa().
.dimsem_unidim <- function(R) {
  k <- ncol(R)
  if (k < 2) {
    return(NA_real_)
  }
  R[!is.finite(R)] <- 0
  off <- upper.tri(R)
  denom <- sum(R[off]^2)
  if (denom < .Machine$double.eps) {
    return(0)
  }

  # Congeneric component: single-factor (minres) implied correlations.
  f <- tryCatch(
    {
      .dimsem_require_namespace("psych")
      fa1 <- suppressWarnings(psych::fa(R, nfactors = 1, fm = "minres",
                                        warnings = FALSE))
      as.numeric(fa1$loadings)
    },
    error = function(e) {
      # Fallback: first eigenvector scaling.
      e1 <- eigen(R, symmetric = TRUE)
      as.numeric(e1$vectors[, 1] * sqrt(max(e1$values[1], 0)))
    }
  )
  Rc <- tcrossprod(f)
  pc <- 1 - sum((R - Rc)[off]^2) / denom

  # Tau-equivalence component: equal loadings, lambda^2 = mean off-diagonal r
  # (floored at 0; a negative manifold cannot be tau-equivalent).
  lam2 <- max(mean(R[off]), 0)
  Rt <- matrix(lam2, k, k)
  tau <- 1 - sum((R - Rt)[off]^2) / denom

  max(min(pc, 1), 0) * max(min(tau, 1), 0)
}

# --- Louvain with "most common" consensus ------------------------------------------

# Louvain with most-common consensus clustering (Lancichinetti & Fortunato,
# 2012; Golino & Christensen, 2025): the algorithm is applied `reps` times
# and the modal solution (up to label permutation) is returned.
#
# `first_pass = TRUE` uses the lower-order (first-pass) memberships of the
# multilevel algorithm (Jimenez et al., 2023), appropriate for LARGE node
# sets where fine-grained communities are wanted. For the small meta-graphs
# dimsem operates on at the hyper level (d <= ~10 nodes), the first pass can
# terminate at dyads and shred genuine communities, so the FINAL modularity
# solution (Blondel et al., 2008) is the default here.
.dimsem_louvain_consensus <- function(g, reps = 2000, resolution = 1,
                                      first_pass = FALSE, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  canon <- function(m) {
    as.integer(factor(m, levels = unique(m)))
  }

  runs <- replicate(reps, {
    cl <- igraph::cluster_louvain(g, weights = igraph::E(g)$weight,
                                  resolution = resolution)
    mm <- tryCatch(cl$memberships, error = function(e) NULL)
    m <- if (isTRUE(first_pass) && !is.null(mm) && is.matrix(mm) &&
             nrow(mm) >= 1) {
      mm[1, ]
    } else {
      as.integer(igraph::membership(cl))
    }
    canon(m)
  })
  runs <- matrix(runs, ncol = reps)

  sigs <- apply(runs, 2, paste, collapse = "-")
  modal_sig <- names(sort(table(sigs), decreasing = TRUE))[1]
  membership <- as.integer(strsplit(modal_sig, "-", fixed = TRUE)[[1]])
  names(membership) <- igraph::V(g)$name

  list(
    membership = membership,
    consensus_share = max(table(sigs)) / reps
  )
}

# --- level-2 network ---------------------------------------------------------------

# EBICglasso on the community-score correlation matrix. An EMPTY estimated
# network is a substantive result (no conditional dependence between
# communities), not an error. Estimation *failure* (small-d instability)
# falls back to lightly ridge-shrunken partial correlations with a warning.
.dimsem_level2_network <- function(scores, glasso_args = list()) {
  .dimsem_require_namespace("qgraph")
  R <- stats::cor(scores, use = "pairwise.complete.obs")
  R[!is.finite(R)] <- 0
  diag(R) <- 1
  n <- nrow(scores)
  args <- .dimsem_merge_args(list(gamma = 0.5), glasso_args)

  net <- tryCatch(
    suppressWarnings(suppressMessages(
      do.call(qgraph::EBICglasso, c(list(S = R, n = n), args))
    )),
    error = function(e) e
  )

  if (inherits(net, "error")) {
    warning("Level-2 EBICglasso failed (", conditionMessage(net),
            "); falling back to ridge-shrunken partial correlations.",
            call. = FALSE)
    Rr <- R + diag(0.05, ncol(R))
    P <- -stats::cov2cor(solve(Rr))
    diag(P) <- 0
    net <- P
  }

  dimnames(net) <- dimnames(R)
  list(network = net, cor = R)
}

# --- strategy 1: TGA-style two-level detection ---------------------------------------

.dimsem_hyper_tga <- function(data, network_matrix, partition,
                              hyper_args = list(), seed = NULL) {
  part <- partition[!is.na(partition)]
  comms <- sort(unique(part))
  d <- length(comms)
  labels <- paste0("F", comms)

  u_hi <- hyper_args$u_hi %||% 0.90
  u_lo <- hyper_args$u_lo %||% 0.50
  reps <- hyper_args$consensus_reps %||% 2000
  resolution <- hyper_args$resolution %||% 1
  min_block_size <- hyper_args$min_block_size %||% 3L

  out <- list(method = "tga", d = d, labels = labels, blocks = list(),
              orphans = character(0), verdict = NA_character_,
              unidim = NA_real_, block_unidim = NULL, level2 = NULL,
              consensus_share = NA_real_, notes = character(0))

  if (d < 2) {
    out$verdict <- "not_applicable"
    out$notes <- "Fewer than two level-1 dimensions; no hyper level exists."
    return(out)
  }

  scores <- .dimsem_network_scores(data, network_matrix, partition)
  lvl2 <- .dimsem_level2_network(scores, hyper_args$glasso_args %||% list())
  out$level2 <- lvl2
  out$unidim <- .dimsem_unidim(lvl2$cor)

  W2 <- abs(lvl2$network)
  diag(W2) <- 0

  if (all(W2 < .Machine$double.eps)) {
    out$verdict <- "none"
    out$notes <- paste0("The level-2 network is empty: community scores show ",
                        "no conditional dependence (unidim u = ",
                        round(out$unidim, 3), ").")
    return(out)
  }

  g2 <- igraph::graph_from_adjacency_matrix(W2, mode = "undirected",
                                            weighted = TRUE, diag = FALSE)
  isolated <- labels[igraph::degree(g2) == 0]
  cons <- .dimsem_louvain_consensus(
    g2, reps = reps, resolution = resolution,
    first_pass = hyper_args$louvain_first_pass %||% FALSE, seed = seed)
  out$consensus_share <- cons$consensus_share

  membership <- cons$membership
  membership[isolated] <- NA_integer_

  blocks <- split(names(membership)[!is.na(membership)],
                  membership[!is.na(membership)])
  block_sizes <- vapply(blocks, length, integer(1))
  out$orphans <- c(isolated, unlist(blocks[block_sizes < 2]))
  blocks <- blocks[block_sizes >= 2]

  # Per-block unidimensionality check; blocks failing clearly are dissolved.
  keep <- rep(TRUE, length(blocks))
  bu <- numeric(length(blocks))
  for (b in seq_along(blocks)) {
    bu[b] <- .dimsem_unidim(lvl2$cor[blocks[[b]], blocks[[b]], drop = FALSE])
    if (is.finite(bu[b]) && bu[b] <= u_lo) {
      keep[b] <- FALSE
    }
  }
  if (any(!keep)) {
    out$orphans <- c(out$orphans, unlist(blocks[!keep]))
    out$notes <- c(out$notes, paste0(
      "Dissolved candidate block(s) with unidim u <= ", u_lo, ": ",
      paste(vapply(blocks[!keep], paste, character(1), collapse = "+"),
            collapse = "; "), "."))
  }
  out$block_unidim <- stats::setNames(bu[keep], vapply(blocks[keep], paste,
                                                       character(1), collapse = "+"))
  blocks <- blocks[keep]

  # Identification gate (see the permutation strategy for the rationale):
  # hyper-factors need at least `min_block_size` (default 3) level-1
  # dimensions; smaller validated groupings are reported as correlated
  # level-1 dimensions rather than hyper-factors.
  small <- vapply(blocks, length, integer(1)) < min_block_size
  if (any(small)) {
    out$notes <- c(out$notes, paste0(
      "Level-2 grouping(s) below min_block_size = ", min_block_size,
      " reported as correlated level-1 dimensions, not hyper-factors: ",
      paste(vapply(blocks[small], paste, character(1), collapse = "+"),
            collapse = "; "), "."))
    out$orphans <- c(out$orphans, unlist(blocks[small]))
    blocks <- blocks[!small]
  }
  # NB: paste0("G", integer(0)) returns "G" (length 1) because paste0()
  # treats zero-length arguments as "" (recycle0 = FALSE); guard the
  # empty case explicitly.
  if (length(blocks) > 0) {
    names(blocks) <- paste0("G", seq_along(blocks))
  }
  out$blocks <- blocks

  n_blocks <- length(blocks)
  if (n_blocks == 0) {
    out$verdict <- "none"
  } else if (n_blocks == 1 && length(blocks[[1]]) == d) {
    # Single community spanning everything: the unidim gate decides whether
    # a single top-level dimension is statistically supported (Samo et al.,
    # 2025, Step 7).
    out$verdict <- if (out$unidim >= u_hi) {
      "single"
    } else if (out$unidim <= u_lo) {
      "none"
    } else {
      "inconclusive"
    }
    if (identical(out$verdict, "none")) {
      out$blocks <- list()
      out$orphans <- labels
      out$notes <- c(out$notes, paste0(
        "Louvain merges all dimensions, but unidim (u = ",
        round(out$unidim, 3), " <= ", u_lo,
        ") rejects a single top-level dimension."))
    }
    if (identical(out$verdict, "inconclusive")) {
      out$notes <- c(out$notes, paste0(
        "unidim u = ", round(out$unidim, 3), " lies in the gray zone (",
        u_lo, ", ", u_hi, "); inspect the level-2 network before imposing ",
        "a hyper-factor."))
    }
  } else {
    out$verdict <- if (n_blocks == 1) "single" else "multiple"
  }

  out
}

# --- strategy 2: permutation-validated higher-order Louvain ---------------------------

.dimsem_hyper_permutation <- function(data, partition,
                                      hyper_args = list(), seed = NULL,
                                      parallel = FALSE, ncores = NULL,
                                      verbose = FALSE) {
  .dimsem_require_namespace("qgraph")

  part <- partition[!is.na(partition)]
  comms <- sort(unique(part))
  d <- length(comms)
  labels <- paste0("F", comms)

  alpha <- hyper_args$alpha %||% 0.05
  # The minimum attainable permutation p-value is 1/(n_perm + 1), so after
  # BH correction the best possible adjusted p for a SINGLE validated pair
  # is n_pairs / (n_perm + 1). The default number of permutations therefore
  # scales with the number of community pairs so that even one lone true
  # connection can survive FDR at `alpha` (factor 2 for headroom),
  # floored at 2000.
  n_pairs <- d * (d - 1) / 2
  n_perm <- hyper_args$n_perm %||% max(2000L, ceiling(2 * n_pairs / alpha))
  reps <- hyper_args$consensus_reps %||% 2000
  resolution <- hyper_args$resolution %||% 1
  # Identification / falsifiability rule: a hyper-factor over two level-1
  # dimensions is an exact reparameterization of their correlation (0 df,
  # untestable); three children are the smallest refutable hyper-dimension
  # (latent-level analogue of the three-indicator rule; cf. Gorsuch, 1997,
  # and the two-node-community handling in Samo et al., 2025).
  min_block_size <- hyper_args$min_block_size %||% 3L
  # Family-resemblance membership rule: a member must have permutation-
  # validated MARGINAL dependence with "majority" (default) or "all" of its
  # co-members. Marginal (not conditional) statistics are used here because
  # EBICglasso legitimately prunes within-family edges via indirect paths,
  # while family resemblance is a claim about raw pairwise relations.
  member_rule <- match.arg(hyper_args$member_rule %||% "majority",
                           choices = c("majority", "all"))
  glasso_args <- .dimsem_merge_args(list(gamma = 0.5),
                                    hyper_args$glasso_args %||% list())

  out <- list(method = "permutation", d = d, labels = labels,
              blocks = list(), orphans = character(0),
              verdict = NA_character_, global_p = NA_real_,
              pair_table = NULL, validated_pairs = NULL,
              evicted = character(0), consensus_share = NA_real_,
              n_perm = n_perm, min_block_size = min_block_size,
              member_rule = member_rule, notes = character(0))

  if (d < 2) {
    out$verdict <- "not_applicable"
    out$notes <- "Fewer than two level-1 dimensions; no hyper level exists."
    return(out)
  }

  items_by_comm <- lapply(comms, function(cc) names(part)[part == cc])
  X <- as.matrix(data[, names(part), drop = FALSE])
  storage.mode(X) <- "numeric"
  n <- nrow(X)

  # Two test statistics per community pair, from the SAME pipeline for
  # observed data and every permutation:
  #   S_cd (conditional): between-community edge mass of the EBICglasso
  #                       network -- drives block DISCOVERY (Stage A).
  #   M_cd (marginal)   : between-community sum of |r| -- drives MEMBERSHIP
  #                       validation (Stage B).
  pair_stats <- function(M) {
    R <- stats::cor(M, use = "pairwise.complete.obs")
    R[!is.finite(R)] <- 0
    diag(R) <- 1
    net <- tryCatch(
      suppressWarnings(suppressMessages(
        do.call(qgraph::EBICglasso, c(list(S = R, n = n), glasso_args))
      )),
      error = function(e) NULL
    )
    if (is.null(net)) {
      net <- R  # degenerate fallback; identical pipeline across permutations
    }
    S <- matrix(0, d, d)
    Mg <- matrix(0, d, d)
    for (a in seq_len(d - 1)) {
      for (b in (a + 1):d) {
        S[a, b] <- sum(abs(net[items_by_comm[[a]], items_by_comm[[b]]]))
        Mg[a, b] <- sum(abs(R[items_by_comm[[a]], items_by_comm[[b]]]))
      }
    }
    list(S = S, M = Mg)
  }

  obs <- pair_stats(X)
  S_obs <- obs$S
  M_obs <- obs$M

  G_obs <- sum(S_obs)

  # --- Permutation index precomputation -----------------------------------
  # RNG REPLICATION GUARANTEE: all permutation indices are drawn SERIALLY
  # from the master RNG stream, in exactly the order the original serial
  # loop drew them (for b in 1:n_perm, then for cc in 1:d, one
  # sample.int(n) each). Workers therefore perform NO random number
  # generation whatsoever -- they are pure functions of the indices handed
  # to them. Results are consequently bit-identical to the serial
  # implementation for a given `seed`, irrespective of `parallel`,
  # `ncores`, scheduling order, or backend.
  if (!is.null(seed)) {
    set.seed(seed)
  }
  perm_idx <- lapply(seq_len(n_perm), function(b) {
    lapply(seq_len(d), function(cc) sample.int(n))
  })

  # One permutation replicate: rebuild the shuffled matrix from stored
  # indices and compute the statistics. Identical arithmetic to the serial
  # loop body.
  perm_one <- function(idx) {
    Xp <- X
    for (cc in seq_len(d)) {
      # Independent row shuffles per community block: preserves each
      # community's internal joint distribution exactly, destroys all
      # between-community dependence (the H0 of mutual independence).
      Xp[, items_by_comm[[cc]]] <- X[idx[[cc]], items_by_comm[[cc]],
                                     drop = FALSE]
    }
    null_b <- pair_stats(Xp)
    list(S = null_b$S >= S_obs,
         M = null_b$M >= M_obs,
         G = sum(null_b$S) >= G_obs)
  }

  use_par <- isTRUE(parallel) &&
    requireNamespace("foreach", quietly = TRUE) &&
    requireNamespace("doSNOW", quietly = TRUE) &&
    requireNamespace("parallel", quietly = TRUE)

  if (isTRUE(parallel) && !use_par) {
    warning("`parallel_hyper = TRUE` requires the `foreach`, `doSNOW`, and ",
            "`parallel` packages; falling back to serial evaluation.",
            call. = FALSE)
  }

  if (use_par) {
    ncores <- .dimsem_resolve_ncores(ncores, n_tasks = n_perm)
  }
  use_par <- use_par && ncores > 1L

  if (use_par) {
    if (isTRUE(verbose)) {
      message("Permutation test: ", n_perm, " replicates on ", ncores,
              " cores (doSNOW).")
    }
    cl <- parallel::makeCluster(ncores)
    on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
    doSNOW::registerDoSNOW(cl)
    # Workers need only qgraph loaded; all other objects travel via the
    # foreach export mechanism. No RNG state is used or set on workers.
    parallel::clusterCall(cl, function() {
      suppressMessages(requireNamespace("qgraph", quietly = TRUE))
      NULL
    })
    b <- NULL  # silence R CMD check note for the foreach iterator
    res_list <- foreach::foreach(
      b = seq_len(n_perm),
      .packages = "qgraph",
      .export = character(0),
      .noexport = character(0)
    ) %dopar% {
      perm_one(perm_idx[[b]])
    }
    # Restore the sequential backend so a later call is unaffected.
    parallel::stopCluster(cl)
    on.exit()
    if (requireNamespace("foreach", quietly = TRUE)) {
      foreach::registerDoSEQ()
    }
  } else {
    if (isTRUE(verbose)) {
      message("Permutation test: ", n_perm, " replicates (serial).")
    }
    res_list <- lapply(perm_idx, perm_one)
  }

  # Reduction is order-independent (sums of indicator matrices), so the
  # accumulated counts do not depend on the order in which replicates
  # complete.
  exceed_S <- matrix(0, d, d)
  exceed_M <- matrix(0, d, d)
  exceed_global <- 0
  for (r in res_list) {
    exceed_S <- exceed_S + r$S
    exceed_M <- exceed_M + r$M
    exceed_global <- exceed_global + r$G
  }

  p_S <- (1 + exceed_S) / (1 + n_perm)
  p_M <- (1 + exceed_M) / (1 + n_perm)
  out$global_p <- (1 + exceed_global) / (1 + n_perm)

  ut <- which(upper.tri(p_S), arr.ind = TRUE)
  pair_table <- data.frame(
    from = labels[ut[, 1]],
    to = labels[ut[, 2]],
    S = S_obs[ut],
    p_S = p_S[ut],
    M = M_obs[ut],
    p_M = p_M[ut],
    stringsAsFactors = FALSE
  )
  pair_table$p_S_fdr <- stats::p.adjust(pair_table$p_S, method = "BH")
  pair_table$p_M_fdr <- stats::p.adjust(pair_table$p_M, method = "BH")
  pair_table$validated <- pair_table$p_S_fdr < alpha
  pair_table$validated_marginal <- pair_table$p_M_fdr < alpha
  out$pair_table <- pair_table[order(pair_table$p_S_fdr, -pair_table$S), ]

  # Symmetric lookup matrices for the membership stage.
  vM <- matrix(FALSE, d, d, dimnames = list(labels, labels))
  mM <- matrix(0, d, d, dimnames = list(labels, labels))
  for (r in seq_len(nrow(pair_table))) {
    i <- pair_table$from[r]; j <- pair_table$to[r]
    vM[i, j] <- vM[j, i] <- pair_table$validated_marginal[r]
    mM[i, j] <- mM[j, i] <- pair_table$M[r]
  }

  if (out$global_p >= alpha) {
    out$verdict <- "none"
    out$orphans <- labels
    out$notes <- paste0(
      "No between-community dependence beyond the mutual-independence ",
      "null (global permutation p = ", round(out$global_p, 4), ").")
    return(out)
  }
  if (!any(pair_table$validated)) {
    out$verdict <- "inconclusive"
    out$orphans <- labels
    out$notes <- paste0(
      "Global between-community dependence is significant (p = ",
      round(out$global_p, 4), "), but no individual pair survives ",
      "BH-FDR at alpha = ", alpha, " with n_perm = ", n_perm,
      "; the dependence is too diffuse to localize. Consider increasing ",
      "hyper_args$n_perm.")
    return(out)
  }

  # --- Stage A: block discovery on the validated CONDITIONAL meta-graph ---
  A <- matrix(0, d, d, dimnames = list(labels, labels))
  val <- pair_table[pair_table$validated, , drop = FALSE]
  for (r in seq_len(nrow(val))) {
    A[val$from[r], val$to[r]] <- val$S[r]
    A[val$to[r], val$from[r]] <- val$S[r]
  }
  g_meta <- igraph::graph_from_adjacency_matrix(A, mode = "undirected",
                                                weighted = TRUE, diag = FALSE)
  isolated <- labels[igraph::degree(g_meta) == 0]

  connected <- setdiff(labels, isolated)
  if (length(connected) <= 2) {
    blocks <- if (length(connected) == 2) list(G1 = connected) else list()
    out$consensus_share <- 1
  } else {
    cons <- .dimsem_louvain_consensus(
      g_meta, reps = reps, resolution = resolution,
      first_pass = hyper_args$louvain_first_pass %||% FALSE, seed = seed)
    out$consensus_share <- cons$consensus_share
    membership <- cons$membership
    membership[isolated] <- NA_integer_
    blocks <- split(names(membership)[!is.na(membership)],
                    membership[!is.na(membership)])
    sizes <- vapply(blocks, length, integer(1))
    isolated <- c(isolated, unlist(blocks[sizes < 2]))
    blocks <- blocks[sizes >= 2]
    if (length(blocks) > 0) {  # paste0("G", integer(0)) returns "G"; guard
      names(blocks) <- paste0("G", seq_along(blocks))
    }
  }

  # --- Stage B: family-resemblance membership validation ------------------
  # Backward elimination: a member needs permutation-validated MARGINAL
  # dependence with `member_rule` of its co-members; the weakest-attached
  # member (fewest validated links, ties broken by smallest marginal mass)
  # is evicted first, and validation is re-assessed on the reduced block.
  # Elimination is sequential/greedy (a purification device built on exact
  # per-step tests), and only ever REMOVES members, which is conservative
  # with respect to hyper-factor imposition.
  evicted <- character(0)
  pruned_blocks <- list()
  for (bl in blocks) {
    members <- bl
    repeat {
      k <- length(members)
      if (k < 2) break
      links <- vapply(members, function(s) {
        sum(vM[s, setdiff(members, s)])
      }, numeric(1))
      need <- if (identical(member_rule, "all")) {
        k - 1
      } else {
        floor((k - 1) / 2) + 1
      }
      ok <- links >= need
      if (all(ok)) break
      mass <- vapply(members, function(s) {
        sum(mM[s, setdiff(members, s)])
      }, numeric(1))
      worst <- members[order(links, mass)][1]
      evicted <- c(evicted, worst)
      members <- setdiff(members, worst)
    }
    if (length(members) >= 2) {
      pruned_blocks <- c(pruned_blocks, list(members))
    } else {
      evicted <- c(evicted, members)
    }
  }
  if (length(evicted) > 0) {
    out$evicted <- evicted
    out$notes <- c(out$notes, paste0(
      "Stray dimension(s) evicted from candidate hyper-communities for ",
      "lacking ", member_rule, "-validated marginal dependence with their ",
      "co-members: ", paste(evicted, collapse = ", "), "."))
  }

  # --- Stage C: identification gate (min_block_size) ----------------------
  final_blocks <- list()
  dyads <- list()
  for (bl in pruned_blocks) {
    if (length(bl) >= min_block_size) {
      final_blocks <- c(final_blocks, list(bl))
    } else {
      dyads <- c(dyads, list(bl))
    }
  }
  if (length(final_blocks) > 0) {  # paste0 zero-length guard, as above
    names(final_blocks) <- paste0("G", seq_along(final_blocks))
  }

  in_block <- unlist(final_blocks)
  out$blocks <- final_blocks
  out$orphans <- setdiff(labels, in_block)

  # Validated dependencies that do NOT justify a hyper-factor are reported,
  # not discarded: a two-dimension "hyper-factor" is an unfalsifiable
  # reparameterization of the pair correlation, so such pairs should be
  # modeled as correlated level-1 dimensions instead.
  vp <- pair_table[pair_table$validated, c("from", "to", "S", "M"),
                   drop = FALSE]
  vp <- vp[!(vp$from %in% in_block & vp$to %in% in_block), , drop = FALSE]
  if (nrow(vp) > 0) {
    out$validated_pairs <- vp
    out$notes <- c(out$notes, paste0(
      "Validated between-dimension dependence outside any identified ",
      "hyper-community (model these as correlated level-1 dimensions, not ",
      "hyper-factors): ",
      paste(paste0(vp$from, "--", vp$to), collapse = ", "), "."))
  }

  out$verdict <- if (length(final_blocks) == 0) {
    "none"
  } else if (length(final_blocks) == 1) {
    "single"
  } else {
    "multiple"
  }

  out
}

# --- orchestrator ------------------------------------------------------------------

.dimsem_hyper_proposal <- function(data, network_matrix, partition,
                                   hyper_args = list(), seed = NULL,
                                   verbose = TRUE, parallel = FALSE,
                                   ncores = NULL) {
  method <- hyper_args$method %||% "both"
  method <- match.arg(method, choices = c("both", "tga", "permutation"))

  results <- list()
  if (method %in% c("both", "tga")) {
    if (isTRUE(verbose)) {
      message("Evaluating hyper-structure: TGA-style level-2 analysis...")
    }
    results$tga <- .dimsem_hyper_tga(data, network_matrix, partition,
                                     hyper_args = hyper_args, seed = seed)
  }
  if (method %in% c("both", "permutation")) {
    if (isTRUE(verbose)) {
      message("Evaluating hyper-structure: permutation-validated ",
              "higher-order Louvain...")
    }
    results$permutation <- .dimsem_hyper_permutation(
      data, partition, hyper_args = hyper_args, seed = seed,
      parallel = parallel, ncores = ncores, verbose = verbose)
  }

  fmt_blocks <- function(res) {
    if (length(res$blocks) == 0) {
      return("-")
    }
    paste(vapply(seq_along(res$blocks), function(i) {
      paste0(names(res$blocks)[i], " = {",
             paste(res$blocks[[i]], collapse = ", "), "}")
    }, character(1)), collapse = "; ")
  }

  verdicts <- vapply(results, function(r) r$verdict, character(1))
  agree <- length(unique(verdicts)) == 1

  # "inconclusive" is a soft verdict (gray-zone unidim, or diffuse
  # dependence that fails to localize); when one strategy is decisive and
  # the other merely inconclusive, follow the decisive one with a caveat.
  decisive <- verdicts[verdicts != "inconclusive"]
  effective <- if (agree) {
    verdicts[[1]]
  } else if (length(decisive) > 0 && length(unique(decisive)) == 1) {
    unique(decisive)
  } else {
    NA_character_
  }
  ref <- if (!is.na(effective) && any(verdicts == effective)) {
    results[[which(verdicts == effective)[1]]]
  } else {
    results[[1]]
  }
  caveat <- if (!agree && !is.na(effective)) {
    inc <- names(verdicts)[verdicts == "inconclusive"]
    paste0(" (The ", paste(inc, collapse = ", "), " strategy was ",
           "inconclusive; this recommendation follows the decisive ",
           "criterion.)")
  } else {
    ""
  }

  recommendation <- if (all(verdicts == "not_applicable")) {
    "Only one level-1 dimension was proposed; there is no hyper level."
  } else if (!is.na(effective) && effective == "none") {
    vp <- results$permutation$validated_pairs
    if (!is.null(vp) && nrow(vp) > 0) {
      paste0("No hyper-factor is identified: validated between-dimension ",
             "dependence exists (",
             paste(paste0(vp$from, "--", vp$to), collapse = ", "),
             "), but no grouping reaches the minimum of ",
             results$permutation$min_block_size %||% 3,
             " member dimensions. Model these dependencies as correlated ",
             "level-1 dimensions: a hyper-factor over fewer than three ",
             "children is an unfalsifiable reparameterization of their ",
             "correlation(s).", caveat)
    } else {
      paste0("Neither imposing a single nor multiple hyper-factors is ",
             "supported: the level-1 dimensions behave as mutually ",
             "independent (personality-like) structures.", caveat)
    }
  } else if (!is.na(effective) && effective == "single") {
    paste0("A single hyper-factor over ",
           fmt_blocks(ref),
           if (length(ref$orphans) > 0) {
             paste0(" is supported, with ",
                    paste(ref$orphans, collapse = ", "),
                    " left as stand-alone dimension(s).")
           } else {
             " is supported (family-resemblance / ideology-like structure)."
           }, caveat)
  } else if (!is.na(effective) && effective == "multiple") {
    paste0("Multiple hyper-factors are supported: ",
           fmt_blocks(ref), ".", caveat)
  } else if (!is.na(effective) && effective == "inconclusive") {
    paste0("Both strategies are inconclusive; inspect `$hyper` details ",
           "and consider increasing hyper_args$n_perm or examining the ",
           "level-2 network directly.")
  } else {
    paste0("The two strategies disagree (",
           paste(names(verdicts), verdicts, sep = ": ", collapse = "; "),
           "); inspect `$hyper` details -- the permutation test is the ",
           "stricter inferential criterion, while the TGA result reflects ",
           "the score-level community structure.")
  }

  list(
    method = method,
    tga = results$tga,
    permutation = results$permutation,
    agreement = agree,
    effective_verdict = effective,
    recommendation = recommendation
  )
}


# --- dimsem_estimate helpers -------------------------------------------------------

# Validate a dimsem_proposal against the supplied data and extract the
# estimation-relevant specification. Hyper blocks are taken with
# permutation-strategy priority (the stricter inferential criterion),
# falling back to TGA -- via .dimsem_hyper_plot_result(), which encodes
# exactly that priority.
.dimsem_extract_proposal <- function(proposal, data, covariates = NULL,
                                     include_hyper = TRUE) {
  # Accept both historical spellings of the class attribute.
  if (!inherits(proposal, "DimSEM_proposal") &&
      !inherits(proposal, "DimSEM_proposal")) {
    stop("`proposal` must be a proposal object from DimSEM_proposal().",
         call. = FALSE)
  }
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("`data` must be a data frame or matrix.", call. = FALSE)
  }
  data <- as.data.frame(data)

  item_names <- proposal$data_info$item_names
  missing_items <- setdiff(item_names, colnames(data))
  if (length(missing_items) > 0) {
    stop("`data` is missing analyzed item column(s): ",
         paste(missing_items, collapse = ", "), call. = FALSE)
  }
  if (nrow(data) != proposal$data_info$n) {
    warning("`data` has ", nrow(data), " rows but the proposal was derived ",
            "from ", proposal$data_info$n, " rows; make sure this is the ",
            "same data set.", call. = FALSE)
  }
  if (!is.null(covariates)) {
    missing_cov <- setdiff(covariates, colnames(data))
    if (length(missing_cov) > 0) {
      stop("Covariate column(s) not found in `data`: ",
           paste(missing_cov, collapse = ", "), call. = FALSE)
    }
    non_num <- covariates[!vapply(data[covariates], is.numeric, logical(1))]
    if (length(non_num) > 0) {
      stop("Covariate(s) must be numeric (dummy-code factors first): ",
           paste(non_num, collapse = ", "), call. = FALSE)
    }
  }

  partition <- proposal$selected$partition
  measurement <- proposal$selected$model_syntax
  if (is.null(measurement) || !nzchar(measurement)) {
    stop("The proposal contains an empty measurement model (no items ",
         "assigned to any dimension); nothing to estimate.", call. = FALSE)
  }
  latent_names <- paste0("F", sort(unique(partition[!is.na(partition)])))

  hyper_blocks <- list()
  hyper_source <- NA_character_
  if (isTRUE(include_hyper)) {
    hp <- .dimsem_hyper_plot_result(proposal$hyper)
    if (!is.null(hp) && length(hp$blocks) > 0) {
      # Keep only blocks whose children exist among the latent dimensions.
      hyper_blocks <- lapply(hp$blocks, intersect, y = latent_names)
      hyper_blocks <- hyper_blocks[vapply(hyper_blocks, length,
                                          integer(1)) >= 2]
      hyper_source <- hp$source
    }
  }

  list(
    item_names = item_names,
    partition = partition,
    latent_names = latent_names,
    measurement_syntax = measurement,
    hyper_blocks = hyper_blocks,
    hyper_source = hyper_source,
    source = proposal$selected$source,
    n = nrow(data)
  )
}

# Assemble full lavaan model syntax: measurement + hyper level + structural
# covariate regressions.
.dimsem_assemble_syntax <- function(spec, covariates = NULL,
                                    covariate_targets = "all",
                                    hyper_loadings = "unity") {
  parts <- c("# Measurement model (from dimsem_propose)",
             spec$measurement_syntax)

  hyper_names <- names(spec$hyper_blocks)
  if (length(spec$hyper_blocks) > 0) {
    parts <- c(parts, "", paste0(
      "# Hyper-factor level (",
      if (identical(hyper_loadings, "unity")) {
        "unit loadings, free variance"
      } else {
        "free loadings, variance fixed via std.lv"
      }, ")"))
    for (g in hyper_names) {
      children <- spec$hyper_blocks[[g]]
      mode_g <- hyper_loadings
      if (identical(mode_g, "free") && length(children) < 3) {
        warning("Hyper-factor ", g, " has only ", length(children),
                " sub-dimension(s); free hyper loadings are not ",
                "identified with fewer than three children. Falling back ",
                "to unity loadings for this block.", call. = FALSE)
        mode_g <- "unity"
      }
      if (identical(mode_g, "unity")) {
        parts <- c(parts,
                   paste0(g, " =~ ", paste(paste0("1*", children),
                                           collapse = " + ")),
                   paste0(g, " ~~ NA*", g))
      } else {
        # std.lv = TRUE fixes var(G) = 1 and frees all loadings.
        parts <- c(parts,
                   paste0(g, " =~ ", paste(children, collapse = " + ")))
      }
    }
  }

  if (!is.null(covariates) && length(covariates) > 0) {
    all_latents <- c(spec$latent_names, hyper_names)
    targets <- if (identical(covariate_targets, "all")) {
      stats::setNames(rep(list(covariates), length(all_latents)), all_latents)
    } else {
      if (!is.list(covariate_targets) || is.null(names(covariate_targets))) {
        stop("`covariate_targets` must be \"all\" or a named list mapping ",
             "latent names to covariate vectors.", call. = FALSE)
      }
      unknown <- setdiff(names(covariate_targets), all_latents)
      if (length(unknown) > 0) {
        stop("`covariate_targets` refers to unknown latent(s): ",
             paste(unknown, collapse = ", "), ". Available: ",
             paste(all_latents, collapse = ", "), call. = FALSE)
      }
      covariate_targets
    }
    parts <- c(parts, "", "# Structural model (covariate regressions)")
    for (lv in names(targets)) {
      cs <- intersect(targets[[lv]], covariates)
      if (length(cs) > 0) {
        parts <- c(parts, paste0(lv, " ~ ", paste(cs, collapse = " + ")))
      }
    }
  }

  paste(parts, collapse = "\n")
}

# Classify and (optionally) standardize external covariates. Continuous
# and ordered-categorical covariates (more than two distinct observed
# values) are z-standardized; binary covariates are kept in their raw
# metric so coefficients remain interpretable as group contrasts. The
# classification is returned alongside the data and reused by the Bayesian
# engine to assign per-covariate priors.
.dimsem_scale_covariates <- function(data, covariates, standardize = TRUE) {
  info <- data.frame(
    covariate = covariates,
    n_values = NA_integer_,
    type = NA_character_,
    standardized = FALSE,
    stringsAsFactors = FALSE
  )
  for (i in seq_along(covariates)) {
    v <- data[[covariates[i]]]
    nv <- length(unique(stats::na.omit(v)))
    info$n_values[i] <- nv
    if (nv <= 2L) {
      info$type[i] <- "binary"
    } else {
      info$type[i] <- "continuous/ordered"
      if (isTRUE(standardize)) {
        data[[covariates[i]]] <- as.numeric(scale(v))  # scale() keeps NAs
        info$standardized[i] <- TRUE
      }
    }
  }
  list(data = data, info = info)
}

# Lightweight stage progress bar with elapsed-time reporting.
.dimsem_stage_progress <- function(stages, enabled = TRUE) {
  n <- length(stages)
  i <- 0L
  t0 <- Sys.time()
  bar <- if (enabled) {
    utils::txtProgressBar(min = 0, max = n, style = 3)
  } else {
    NULL
  }
  list(
    tick = function() {
      i <<- i + 1L
      if (!is.null(bar)) {
        utils::setTxtProgressBar(bar, i)
      }
    },
    close = function() {
      if (!is.null(bar)) {
        close(bar)
        message("Done in ",
                format(round(as.numeric(difftime(Sys.time(), t0,
                                                 units = "secs")), 1)),
                " s.")
      }
    }
  )
}

# ML engine: lavaan::sem with any lavaan estimator and missing-data method.
# For larger data sets a quick pilot fit on a subsample provides an
# approximate time-to-completion estimate (lavaan itself exposes no
# iteration callback, so within-fit progress cannot be shown).
.dimsem_estimate_ml <- function(syntax, data, estimator = "MLR",
                                missing = "ml", lavaan_args = list(),
                                progress = TRUE, seed = NULL,
                                verbose = TRUE) {
  .dimsem_require_namespace("lavaan")
  if (!is.null(seed)) {
    set.seed(seed)
  }

  defaults <- list(std.lv = TRUE)
  args <- .dimsem_merge_args(defaults, lavaan_args)
  args$estimator <- args$estimator %||% estimator
  args$missing <- args$missing %||% missing

  ml_family <- toupper(args$estimator) %in%
    c("ML", "MLR", "MLM", "MLMV", "MLMVS", "MLF")
  if (!ml_family && tolower(args$missing) %in%
      c("ml", "fiml", "ml.x", "direct")) {
    warning("FIML missing handling requires an ML-family estimator; lavaan ",
            "will likely switch or error with estimator = \"",
            args$estimator, "\". Consider missing = \"pairwise\".",
            call. = FALSE)
  }

  # Pilot-based ETA: fit on a subsample and extrapolate roughly linearly
  # in n (a heuristic -- optimizer path length varies with the data).
  n <- nrow(data)
  if (isTRUE(progress) && n > 600) {
    n_pilot <- 300L
    pilot_rows <- seq_len(n_pilot)
    t_pilot <- system.time(
      tryCatch(suppressWarnings(.dimsem_lavaan_call(
        syntax, data[pilot_rows, , drop = FALSE], args)),
        error = function(e) NULL)
    )[["elapsed"]]
    if (is.finite(t_pilot) && t_pilot > 0) {
      eta <- t_pilot * (n / n_pilot)
      message("Estimated time to completion: ~",
              format(round(eta, 1)), " s (pilot extrapolation).")
    }
  }

  fit <- tryCatch(.dimsem_lavaan_call(syntax, data, args),
                  error = function(e) e)
  if (inherits(fit, "error")) {
    warning("lavaan estimation failed: ", conditionMessage(fit),
            call. = FALSE)
    return(list(fit = NULL, converged = FALSE, factor_cor = NULL,
                structural = NULL, fit_measures = NULL,
                error = conditionMessage(fit)))
  }

  converged <- isTRUE(tryCatch(lavaan::lavInspect(fit, "converged"),
                               error = function(e) FALSE))
  factor_cor <- if (converged) {
    tryCatch(lavaan::lavInspect(fit, "cor.lv"), error = function(e) NULL)
  } else {
    NULL
  }
  structural <- if (converged) {
    tryCatch({
      pe <- lavaan::parameterEstimates(fit, standardized = TRUE)
      pe <- pe[pe$op == "~", c("lhs", "rhs", "est", "se", "z", "pvalue",
                               "std.all"), drop = FALSE]
      names(pe)[1:2] <- c("latent", "covariate")
      pe
    }, error = function(e) NULL)
  } else {
    NULL
  }
  fit_measures <- if (converged) {
    tryCatch({
      fm <- lavaan::fitMeasures(fit)
      plain <- c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr")
      scaled <- paste0(plain, ".scaled")
      picked <- ifelse(scaled %in% names(fm) & !is.na(fm[scaled]),
                       scaled, plain)
      keep <- fm[picked]
      names(keep) <- plain
      attr(keep, "scaled") <- any(picked != plain)
      keep
    }, error = function(e) NULL)
  } else {
    NULL
  }
  if (!converged) {
    warning("The lavaan model did not converge; inspect `$fit`.",
            call. = FALSE)
  }

  list(fit = fit, converged = converged, factor_cor = factor_cor,
       structural = structural, fit_measures = fit_measures, error = NULL)
}

# lavaan::sem via a quoted-symbol call: splicing the data.frame VALUE into
# the call makes lavaan's non-standard evaluation resolve `data` to the
# utils::data() closure (see .dimsem_fit_cfa for the original diagnosis).
.dimsem_lavaan_call <- function(syntax, data, args) {
  lav_call <- as.call(c(
    list(quote(lavaan::sem)),
    list(model = syntax, data = quote(data)),
    args
  ))
  eval(lav_call)
}

# --- Bayesian (Stan) backend --------------------------------------------------------

# Translate the dimsem model family into a Stan program. Deliberately NOT a
# general lavaan-to-Stan compiler (that is blavaan's remit; Merkle &
# Rosseel, 2018, whose parameterization conventions the measurement part
# follows -- see the package citation file): this generator covers exactly
# the models dimsem_propose() emits, which is what makes CPU/GPU-oriented
# code generation tractable:
#   * first-order factors from the item partition (free loadings, one
#     sign-anchored loading per factor; standardized latent disturbances),
#   * optional hyper-factors with loadings fixed to 1 and free variance
#     (family-resemblance convention),
#   * covariates regressed on any subset of latents,
#   * likelihood over OBSERVED item entries only (long format), i.e.
#     full-information under MAR with no case deletion,
#   * `normal` or `ordinal` (ordered-logistic) item likelihood,
#   * within-chain parallelism via reduce_sum over respondents
#     (threads_per_chain > 1), and OpenCL-compatible constructs for GPU
#     acceleration via cmdstanr's opencl_ids.
.dimsem_stan_model <- function(likelihood = c("normal", "ordinal"),
                               n_hyper = 0L, n_cov = 0L,
                               hyper_loadings = "unity",
                               priors = list()) {
  likelihood <- match.arg(likelihood)
  has_hyper <- n_hyper > 0L
  has_cov <- n_cov > 0L
  free_hyper <- has_hyper && identical(hyper_loadings, "free")

  # Prior defaults; any entry can be overridden via bayes_args$priors.
  pr <- .dimsem_merge_args(list(
    lambda = "normal(0, 2)",
    nu = "normal(0, 5)",
    sigma = "student_t(3, 0, 2)",
    tau_g = "student_t(3, 0, 2)",
    delta = "normal(0, 2)",
    beta = "std_normal()",
    beta_binary = "normal(0, 2)",
    cut = "normal(0, 4)"
  ), priors)

  functions_block <- paste0(
    "functions {
  real partial_ll(array[] int slice_obs, int start, int end,
                  vector y_obs, array[] int obs_row, array[] int obs_col,
                  matrix eta, vector lambda, vector nu, vector sigma,
                  array[] int item_factor",
    if (identical(likelihood, "ordinal")) ",
                  array[] vector cut, array[] int n_cat" else "", ") {
    real ll = 0;
    for (s in start:end) {
      int r = obs_row[s];
      int p = obs_col[s];
      real mu = nu[p] + lambda[p] * eta[r, item_factor[p]];
",
    if (identical(likelihood, "normal"))
      "      ll += normal_lpdf(y_obs[s] | mu, sigma[p]);"
    else
      "      ll += ordered_logistic_lpmf(to_int(y_obs[s]) | mu / sigma[p],
                                   head(cut[p], n_cat[p] - 1));",
    "
    }
    return ll;
  }
}
")

  data_block <- paste0(
    "data {
  int<lower=1> N;                    // respondents
  int<lower=1> P;                    // items
  int<lower=1> K;                    // first-order factors
  array[P] int<lower=1, upper=K> item_factor;
  array[K] int<lower=1, upper=P> ref_item;   // sign-reference item per factor
  int<lower=0> n_obs;                // observed item entries (long format)
  array[n_obs] int<lower=1, upper=N> obs_row;
  array[n_obs] int<lower=1, upper=P> obs_col;
  vector[n_obs] y_obs;               // observed responses
",
    if (has_hyper)
      "  int<lower=1> G;                    // hyper-factors
  array[K] int<lower=0, upper=G> factor_hyper;   // 0 = stand-alone
" else "",
    if (has_hyper)
      "  array[G] int<lower=1, upper=K> ref_child;   // sign-reference child per hyper
" else "",
    if (has_cov)
      "  int<lower=1> Q;                    // covariates
  matrix[N, Q] X;
  array[Q] int<lower=0, upper=1> x_binary; // raw binary predictor flag
  array[K] int<lower=0, upper=1> reg_f;    // factor receives regressions?
",
    if (has_cov && has_hyper)
      "  array[G] int<lower=0, upper=1> reg_g;
" else "",
    if (identical(likelihood, "ordinal"))
      "  int<lower=2> max_cat;
  array[P] int<lower=2> n_cat;
" else "",
    "  int<lower=1> grainsize;
}
")

  params_block <- paste0(
    "parameters {
  vector[P] lambda_un;                    // unconstrained; sign handled by
                                          // in-iterations relabeling in the
                                          // generated quantities (blavaan
                                          // convention; Merkle et al., 2021)
  vector[P] nu;                           // item intercepts
  vector<lower=0>[P] sigma;               // item residual scales
  matrix[N, K] z_f;                       // std-normal factor innovations
",
    if (has_hyper)
      "  matrix[N, G] z_g;
" else "",
    if (has_hyper && !free_hyper)
      "  vector<lower=0>[G] tau_g;               // hyper-factor SDs (free)
" else "",
    if (free_hyper)
      "  vector[K] delta_un;                     // unconstrained hyper loadings
" else "",
    if (has_cov)
      "  matrix[Q, K] beta_f;                    // covariate -> factor
" else "",
    if (has_cov && has_hyper)
      "  matrix[Q, G] beta_g;                    // covariate -> hyper-factor
" else "",
    if (identical(likelihood, "ordinal"))
      "  array[P] ordered[max_cat - 1] cut;      // per-item cutpoints
" else "",
    "}
")

  transformed_block <- paste0(
    "transformed parameters {
",
    if (free_hyper)
      "  vector[K] delta_l;
" else "",
    if (has_hyper)
      "  matrix[N, G] g_score;
" else "",
    "  matrix[N, K] eta;
",
    if (free_hyper)
      "  for (k in 1:K) {
    delta_l[k] = factor_hyper[k] == 0 ? 0 : delta_un[k];
  }
" else "",
    if (has_hyper) paste0(
      "  for (g in 1:G) {
    g_score[, g] = ", if (has_cov) "X * beta_g[, g] + " else "",
      if (free_hyper) "z_g[, g];   // var(G) = 1 (free loadings)"
      else "tau_g[g] * z_g[, g];",
      "
  }
"),
    "  for (k in 1:K) {
    eta[, k] = ",
    if (has_cov) "(reg_f[k] == 1 ? X * beta_f[, k] : rep_vector(0, N)) + " else "",
    if (has_hyper && !free_hyper)
      "(factor_hyper[k] > 0 ? g_score[, factor_hyper[k]] : rep_vector(0, N)) + "
    else "",
    if (free_hyper)
      "(factor_hyper[k] > 0 ? delta_l[k] * g_score[, factor_hyper[k]] : rep_vector(0, N)) + "
    else "",
    "z_f[, k];   // hyper path per loading mode; std disturbance (var 1)
  }
}
")

  model_block <- paste0(
    paste0("model {
  // Priors (weakly informative defaults; override via bayes_args$priors)
  lambda_un ~ ", pr$lambda, ";
  nu ~ ", pr$nu, ";
  sigma ~ ", pr$sigma, ";
  to_vector(z_f) ~ std_normal();
"),
    if (has_hyper)
      "  to_vector(z_g) ~ std_normal();
" else "",
    if (has_hyper && !free_hyper) paste0(
      "  tau_g ~ ", pr$tau_g, ";
") else "",
    if (free_hyper) paste0(
      "  delta_un ~ ", pr$delta, ";
") else "",
    if (has_cov) paste0(
      "  // Standard-normal priors for standardized covariates; a wider
  // normal(0, 2) for raw binary predictors (one raw unit ~ two
  // predictor SDs for a balanced binary variable).
  for (q in 1:Q) {
    if (x_binary[q] == 1) {
      beta_f[q] ~ ", pr$beta_binary, ";
    } else {
      beta_f[q] ~ ", pr$beta, ";
    }
  }
") else "",
    if (has_cov && has_hyper) paste0(
      "  for (q in 1:Q) {
    if (x_binary[q] == 1) {
      beta_g[q] ~ ", pr$beta_binary, ";
    } else {
      beta_g[q] ~ ", pr$beta, ";
    }
  }
") else "",
    if (identical(likelihood, "ordinal")) paste0(
      "  for (p in 1:P) cut[p] ~ ", pr$cut, ";
") else "",
    "
  // Full-information likelihood over OBSERVED entries only (MAR),
  // parallelized within-chain over entries via reduce_sum.
  target += reduce_sum(partial_ll, obs_row, grainsize,
                       y_obs, obs_row, obs_col, eta, lambda_un, nu, sigma,
                       item_factor",
    if (identical(likelihood, "ordinal")) ", cut, n_cat" else "", ");
}
")

  gq_block <- paste0(
"generated quantities {
  // In-iterations sign relabeling (blavaan convention; Merkle et al.,
  // 2021; evaluated in Chen, Miocevic, & Falk, 2025): the likelihood is
  // reflection-invariant, so raw *_un draws may occupy mirror modes
  // across chains BY DESIGN. Per iteration, each factor is reflected so
  // its reference item loading is positive, with all sign-bearing
  // quantities of that factor reflected jointly; hyper orientations
  // follow their reference child. Report and diagnose these corrected
  // quantities, never the *_un scaffolding.
  vector[K] sign_f;
  vector[P] lambda;
",
if (has_cov)
"  matrix[Q, K] beta_f_cor;
" else "",
if (has_hyper)
"  vector[G] sign_g;
" else "",
if (free_hyper)
"  vector[K] delta;
  vector[K] delta_std;
" else "",
if (has_cov && has_hyper)
"  matrix[Q, G] beta_g_cor;
" else "",
if (has_hyper && !free_hyper)
"  vector[G] rho_g;
" else "",
"  for (k in 1:K) {
    sign_f[k] = lambda_un[ref_item[k]] < 0 ? -1 : 1;
  }
  for (p in 1:P) {
    lambda[p] = sign_f[item_factor[p]] * lambda_un[p];
  }
",
if (has_cov)
"  for (k in 1:K) {
    beta_f_cor[, k] = sign_f[k] * beta_f[, k];
  }
" else "",
if (free_hyper)
"  {
    vector[K] delta_tmp;
    for (k in 1:K) delta_tmp[k] = sign_f[k] * delta_l[k];
    for (g in 1:G) sign_g[g] = delta_tmp[ref_child[g]] < 0 ? -1 : 1;
    for (k in 1:K) {
      delta[k] = factor_hyper[k] == 0 ? 0
                 : sign_g[factor_hyper[k]] * delta_tmp[k];
      delta_std[k] = delta[k] / sqrt(square(delta[k]) + 1);
    }
  }
" else "",
if (has_hyper && !free_hyper)
"  // Unity loadings: per-factor reflection is only a joint block
  // symmetry, so the hyper orientation follows its reference child.
  for (g in 1:G) sign_g[g] = sign_f[ref_child[g]];
  for (g in 1:G) rho_g[g] = square(tau_g[g]) / (square(tau_g[g]) + 1);
" else "",
if (has_cov && has_hyper)
"  for (g in 1:G) {
    beta_g_cor[, g] = sign_g[g] * beta_g[, g];
  }
" else "",
"}
")

  paste0(functions_block, data_block, params_block, transformed_block,
         model_block, gq_block)
}

# Build the Stan data list from the estimation spec, in the long
# observed-entries format (no case deletion, no imputation).
.dimsem_stan_data <- function(spec, data, covariates = NULL,
                              covariate_targets = "all",
                              likelihood = "normal", grainsize = NULL,
                              hyper_loadings = "unity",
                              covariate_info = NULL,
                              sign_ref = NULL) {
  items <- names(spec$partition)[!is.na(spec$partition)]
  part <- spec$partition[items]
  K <- length(unique(part))
  factor_index <- as.integer(factor(part, levels = sort(unique(part))))

  Y <- as.matrix(data[, items, drop = FALSE])
  storage.mode(Y) <- "numeric"
  obs <- which(!is.na(Y), arr.ind = TRUE)

  # Sign-reference items (one per factor): data-driven default = the
  # marginal "hub" item, i.e. the item with the highest mean absolute
  # Pearson correlation with the other members of its community. A strong
  # reference is what makes relabeling reliable (Chen, Miocevic, & Falk,
  # 2025, show poor references degrade every sign-handling approach; cf.
  # the maximin-|draw| posterior rule of Kastner et al., 2017).
  ref_item <- .dimsem_sign_refs(data, items, factor_index,
                                sign_ref = sign_ref)

  sd <- list(
    N = nrow(Y), P = length(items), K = K,
    item_factor = factor_index, ref_item = ref_item,
    n_obs = nrow(obs), obs_row = unname(obs[, 1]),
    obs_col = unname(obs[, 2]), y_obs = unname(Y[obs]),
    grainsize = grainsize %||% max(1L, floor(nrow(obs) / 40L))
  )

  hyper_names <- names(spec$hyper_blocks)
  if (length(hyper_names) > 0) {
    latent_labels <- paste0("F", sort(unique(part)))
    fh <- integer(K)
    for (gi in seq_along(hyper_names)) {
      fh[latent_labels %in% spec$hyper_blocks[[gi]]] <- gi
    }
    sd$G <- length(hyper_names)
    sd$factor_hyper <- fh
    rc <- integer(length(hyper_names))
    for (gi in seq_along(hyper_names)) {
      kids <- which(fh == gi)
      rc[gi] <- if (length(kids) > 0) kids[1] else 1L
      gl <- hyper_names[gi]
      if (!is.null(sign_ref) && gl %in% names(sign_ref)) {
        cand <- match(sign_ref[[gl]], latent_labels)
        if (!is.na(cand) && fh[cand] == gi) rc[gi] <- cand
      }
    }
    sd$ref_child <- rc
  }

  if (!is.null(covariates) && length(covariates) > 0) {
    X <- as.matrix(data[, covariates, drop = FALSE])
    storage.mode(X) <- "numeric"
    if (anyNA(X)) {
      warning("Covariates contain missing values; rows with missing ",
              "covariates contribute their item information but the ",
              "current Stan program requires complete X -- these NAs are ",
              "mean-imputed IN THE COVARIATES ONLY (flagged here ",
              "explicitly; consider multiple imputation upstream).",
              call. = FALSE)
      for (j in seq_len(ncol(X))) {
        X[is.na(X[, j]), j] <- mean(X[, j], na.rm = TRUE)
      }
    }
    latent_labels <- paste0("F", sort(unique(part)))
    reg_f <- as.integer(latent_labels %in% .dimsem_reg_targets(
      covariate_targets, latent_labels, hyper_names))
    sd$Q <- ncol(X)
    sd$X <- X
    sd$x_binary <- if (!is.null(covariate_info)) {
      as.integer(covariate_info$type[match(covariates,
                                           covariate_info$covariate)] ==
                   "binary")
    } else {
      vapply(covariates, function(cv) {
        as.integer(length(unique(stats::na.omit(data[[cv]]))) <= 2L)
      }, integer(1), USE.NAMES = FALSE)
    }
    sd$reg_f <- reg_f
    if (length(hyper_names) > 0) {
      sd$reg_g <- as.integer(hyper_names %in% .dimsem_reg_targets(
        covariate_targets, latent_labels, hyper_names, hyper = TRUE))
    }
  }

  if (identical(likelihood, "ordinal")) {
    n_cat <- vapply(items, function(it) {
      length(unique(stats::na.omit(data[[it]])))
    }, integer(1))
    sd$max_cat <- max(n_cat)
    sd$n_cat <- n_cat
  }

  sd
}

.dimsem_reg_targets <- function(covariate_targets, latent_labels,
                                hyper_names, hyper = FALSE) {
  pool <- if (hyper) hyper_names else latent_labels
  if (identical(covariate_targets, "all")) {
    return(pool)
  }
  intersect(names(covariate_targets), pool)
}

# Bayesian engine driver. backend = "cmdstanr" compiles and samples
# (with parallel_chains, threads_per_chain, and opencl_ids passed
# through); backend = "dry_run" returns the generated program and data
# without requiring a Stan toolchain.
.dimsem_estimate_bayes <- function(spec, syntax, data, covariates = NULL,
                                   covariate_targets = "all",
                                   covariate_info = NULL,
                                   hyper_loadings = "unity",
                                   bayes_args = list(),
                                   force_recompile = FALSE, seed = NULL,
                                   verbose = TRUE) {
  backend <- bayes_args$backend %||% "cmdstanr"
  likelihood <- match.arg(bayes_args$likelihood %||% "normal",
                          c("normal", "ordinal"))

  # Guard: free hyper loadings need >= 3 children per block (mirrors the
  # lavaan-side guard in .dimsem_assemble_syntax()).
  if (identical(hyper_loadings, "free") &&
      any(vapply(spec$hyper_blocks, length, integer(1)) < 3)) {
    warning("Free hyper loadings require at least three sub-dimensions ",
            "per hyper-factor; using unity loadings for the Stan model.",
            call. = FALSE)
    hyper_loadings <- "unity"
  }

  # GPU-friendly program selection: whenever the user supplies OpenCL
  # device ids (or sets bayes_args$gpu = TRUE), the generator emits a fully
  # vectorized likelihood -- single large sampling statements over the long
  # observed-entries format -- instead of the reduce_sum/scalar-loop CPU
  # variant, and the data builder performs the matching index/segment
  # transformations automatically.
  gpu <- !is.null(bayes_args$opencl_ids) || isTRUE(bayes_args$gpu)
  if (gpu && (bayes_args$threads_per_chain %||% 1L) > 1L) {
    warning("The GPU-vectorized likelihood replaces reduce_sum; ",
            "`threads_per_chain` is ignored when `opencl_ids` is set ",
            "(within-statement OpenCL parallelism takes its place; ",
            "between-chain parallelism via `parallel_chains` still ",
            "applies).", call. = FALSE)
  }

  if (gpu) {
    stan_code <- .dimsem_stan_model_gpu(
      likelihood = likelihood,
      n_hyper = length(spec$hyper_blocks),
      n_cov = length(covariates %||% character(0)),
      hyper_loadings = hyper_loadings,
      priors = bayes_args$priors %||% list()
    )
    stan_data <- .dimsem_stan_data(
      spec, data, covariates = covariates,
      covariate_targets = covariate_targets,
      likelihood = likelihood, grainsize = bayes_args$grainsize,
      hyper_loadings = hyper_loadings, covariate_info = covariate_info,
      sign_ref = bayes_args$sign_ref
    )
    stan_data <- .dimsem_stan_data_gpu(stan_data, likelihood = likelihood)
  } else {
    stan_code <- .dimsem_stan_model(
      likelihood = likelihood,
      n_hyper = length(spec$hyper_blocks),
      n_cov = length(covariates %||% character(0)),
      hyper_loadings = hyper_loadings,
      priors = bayes_args$priors %||% list()
    )
    stan_data <- .dimsem_stan_data(
      spec, data, covariates = covariates,
      covariate_targets = covariate_targets,
      likelihood = likelihood, grainsize = bayes_args$grainsize,
      hyper_loadings = hyper_loadings, covariate_info = covariate_info,
      sign_ref = bayes_args$sign_ref
    )
  }
  stan <- list(model_code = stan_code, data = stan_data,
               likelihood = likelihood, hyper_loadings = hyper_loadings,
               variant = if (gpu) "gpu_vectorized" else "cpu_reduce_sum")

  if (identical(backend, "dry_run") ||
      !requireNamespace("cmdstanr", quietly = TRUE)) {
    if (!identical(backend, "dry_run")) {
      warning("Package `cmdstanr` is not available; returning the generated ",
              "Stan program and data (dry run). Install cmdstanr and ",
              "CmdStan to sample.", call. = FALSE)
    }
    return(list(fit = NULL, converged = FALSE, factor_cor = NULL,
                structural = NULL, fit_measures = NULL, stan = stan,
                param_map = .dimsem_stan_param_map(spec, covariates,
                                                   hyper_loadings),
                proposal_partition = spec$partition, error = NULL))
  }

  chains <- bayes_args$chains %||% 4L
  threads <- if (gpu) 1L else bayes_args$threads_per_chain %||% 1L
  opencl_ids <- bayes_args$opencl_ids
  cpp_options <- list()
  if (threads > 1L) cpp_options$stan_threads <- TRUE
  if (!is.null(opencl_ids)) cpp_options$stan_opencl <- TRUE

  stan_file <- cmdstanr::write_stan_file(stan_code)
  if (isTRUE(force_recompile)) {
    .dimsem_clear_stan_exe(stan_file, verbose = verbose)
  }
  compile_args <- list(stan_file, cpp_options = cpp_options)
  # Recent cmdstanr versions expose force_recompile directly; pass it
  # through where available (the manual deletion above already guarantees
  # the rebuild on older versions).
  if ("force_recompile" %in% names(formals(cmdstanr::cmdstan_model))) {
    compile_args$force_recompile <- force_recompile
  }
  mod <- do.call(cmdstanr::cmdstan_model, compile_args)

  sample_args <- list(
    data = stan_data,
    chains = chains,
    parallel_chains = bayes_args$parallel_chains %||% chains,
    iter_warmup = bayes_args$iter_warmup %||% 1000L,
    iter_sampling = bayes_args$iter_sampling %||% 1000L,
    adapt_delta = bayes_args$adapt_delta %||% 0.9,
    refresh = bayes_args$refresh %||% if (isTRUE(verbose)) 100L else 0L,
    seed = seed %||% sample.int(.Machine$integer.max, 1L)
  )
  if (threads > 1L) sample_args$threads_per_chain <- threads
  if (!is.null(opencl_ids)) sample_args$opencl_ids <- opencl_ids

  fit <- do.call(mod$sample, sample_args)

  diag <- tryCatch(fit$diagnostic_summary(), error = function(e) NULL)
  converged <- !is.null(diag) &&
    sum(diag$num_divergent %||% 0) == 0 &&
    sum(diag$num_max_treedepth %||% 0) == 0

  list(fit = fit, converged = converged, factor_cor = NULL,
       structural = tryCatch(
         as.data.frame(fit$summary(variables = c(
           if (!is.null(covariates)) "beta_f_cor",
           if (!is.null(covariates) && length(spec$hyper_blocks) > 0)
             "beta_g_cor"))),
         error = function(e) NULL),
       fit_measures = NULL, stan = stan,
       param_map = .dimsem_stan_param_map(spec, covariates, hyper_loadings),
       proposal_partition = spec$partition, error = NULL)
}


# GPU-friendly Stan program: fully vectorized likelihood, no reduce_sum.
#
# Rationale (cf. the matrix-based _lpmf pattern: Stan's OpenCL backend accelerates
# LARGE vectorized lpdf/lpmf calls; per-entry scalar calls inside loops (or inside reduce_sum slices)
# trigger no OpenCL kernels and pay CPU<->GPU transfer for nothing. This
# variant therefore expresses the entire observed-entries likelihood as
#   normal : ONE vectorized normal_lpdf over all n_obs entries, with the
#            parameter-dependent mean assembled from vectorized gathers
#            (eta_obs = to_vector(eta)[eta_idx], a single indexed read);
#   ordinal: P vectorized ordered_logistic_lpmf calls, one per item over
#            that item's contiguous (pre-sorted) observation segment, so
#            each call shares a single cutpoint vector.
# The index and segment structures are precomputed in R by
# .dimsem_stan_data_gpu(), keeping the Stan program free of scalar loops
# over observations.
.dimsem_stan_model_gpu <- function(likelihood = c("normal", "ordinal"),
                                   n_hyper = 0L, n_cov = 0L,
                                   hyper_loadings = "unity",
                                   priors = list()) {
  likelihood <- match.arg(likelihood)
  has_hyper <- n_hyper > 0L
  has_cov <- n_cov > 0L
  free_hyper <- has_hyper && identical(hyper_loadings, "free")

  pr <- .dimsem_merge_args(list(
    lambda = "normal(0, 2)",
    nu = "normal(0, 5)",
    sigma = "student_t(3, 0, 2)",
    tau_g = "student_t(3, 0, 2)",
    delta = "normal(0, 2)",
    beta = "std_normal()",
    beta_binary = "normal(0, 2)",
    cut = "normal(0, 4)"
  ), priors)

  data_block <- paste0(
"data {
  int<lower=1> N;
  int<lower=1> P;
  int<lower=1> K;
  array[P] int<lower=1, upper=K> item_factor;
  array[K] int<lower=1, upper=P> ref_item;   // sign-reference item per factor
  int<lower=0> n_obs;
  array[n_obs] int<lower=1, upper=N> obs_row;
  array[n_obs] int<lower=1, upper=P> obs_col;
  array[n_obs] int<lower=1> eta_idx;      // (factor-1)*N + row: gather index
",
if (identical(likelihood, "normal"))
"  vector[n_obs] y_obs;
" else
"  array[n_obs] int y_int;                 // entries sorted by item
  array[P] int<lower=0> seg_start;        // per-item segment bounds
  array[P] int<lower=0> seg_end;          // seg_end < seg_start = no obs
  int<lower=2> max_cat;
  array[P] int<lower=2> n_cat;
",
if (has_hyper)
"  int<lower=1> G;
  array[K] int<lower=0, upper=G> factor_hyper;
" else "",
if (has_hyper)
"  array[G] int<lower=1, upper=K> ref_child;   // sign-reference child per hyper
" else "",
if (has_cov)
"  int<lower=1> Q;
  matrix[N, Q] X;
  array[Q] int<lower=0, upper=1> x_binary;
  array[K] int<lower=0, upper=1> reg_f;
" else "",
if (has_cov && has_hyper)
"  array[G] int<lower=0, upper=1> reg_g;
" else "",
"}
")

  params_block <- paste0(
"parameters {
  vector[P] lambda_un;                    // unconstrained; sign handled by
                                          // in-iterations relabeling in the
                                          // generated quantities (blavaan
                                          // convention; Merkle et al., 2021)
  vector[P] nu;
  vector<lower=0>[P] sigma;
  matrix[N, K] z_f;
",
if (has_hyper)
"  matrix[N, G] z_g;
" else "",
if (has_hyper && !free_hyper)
"  vector<lower=0>[G] tau_g;
" else "",
if (free_hyper)
"  vector[K] delta_un;                     // unconstrained hyper loadings
" else "",
if (has_cov)
"  matrix[Q, K] beta_f;
" else "",
if (has_cov && has_hyper)
"  matrix[Q, G] beta_g;
" else "",
if (identical(likelihood, "ordinal"))
"  array[P] ordered[max_cat - 1] cut;
" else "",
"}
")

  transformed_block <- paste0(
"transformed parameters {
",
if (free_hyper)
"  vector[K] delta_l;
" else "",
if (has_hyper)
"  matrix[N, G] g_score;
" else "",
"  matrix[N, K] eta;
",
if (free_hyper)
"  for (k in 1:K) {
    delta_l[k] = factor_hyper[k] == 0 ? 0 : delta_un[k];
  }
" else "",
if (has_hyper) paste0(
"  for (g in 1:G) {
    g_score[, g] = ", if (has_cov) "X * beta_g[, g] + " else "",
                     if (free_hyper) "z_g[, g];"
                     else "tau_g[g] * z_g[, g];",
"
  }
"),
"  for (k in 1:K) {
    eta[, k] = ",
if (has_cov) "(reg_f[k] == 1 ? X * beta_f[, k] : rep_vector(0, N)) + " else "",
if (has_hyper && !free_hyper)
"(factor_hyper[k] > 0 ? g_score[, factor_hyper[k]] : rep_vector(0, N)) + "
else "",
if (free_hyper)
"(factor_hyper[k] > 0 ? delta_l[k] * g_score[, factor_hyper[k]] : rep_vector(0, N)) + "
else "",
"z_f[, k];
  }
}
")

  likelihood_block <- if (identical(likelihood, "normal")) paste0(
"
  // Fully vectorized likelihood: one OpenCL-eligible statement over all
  // observed entries; the mean is assembled from vectorized gathers only.
  {
    vector[n_obs] eta_obs = to_vector(eta)[eta_idx];
    vector[n_obs] mu = nu[obs_col] + lambda_un[obs_col] .* eta_obs;
    target += normal_lpdf(y_obs | mu, sigma[obs_col]);
  }
") else paste0(
"
  // Per-item vectorized ordered-logistic calls over contiguous pre-sorted
  // segments: each call shares one cutpoint vector, no scalar loops over
  // observations.
  {
    vector[n_obs] eta_obs = to_vector(eta)[eta_idx];
    for (p in 1:P) {
      if (seg_end[p] >= seg_start[p]) {
        int a = seg_start[p];
        int b = seg_end[p];
        target += ordered_logistic_lpmf(
          y_int[a:b] |
          (nu[p] + lambda_un[p] * eta_obs[a:b]) / sigma[p],
          head(cut[p], n_cat[p] - 1));
      }
    }
  }
")

  model_block <- paste0(
"model {
  lambda_un ~ ", pr$lambda, ";
  nu ~ ", pr$nu, ";
  sigma ~ ", pr$sigma, ";
  to_vector(z_f) ~ std_normal();
",
if (has_hyper)
"  to_vector(z_g) ~ std_normal();
" else "",
if (has_hyper && !free_hyper) paste0("  tau_g ~ ", pr$tau_g, ";
") else "",
if (free_hyper) paste0("  delta_un ~ ", pr$delta, ";
") else "",
if (has_cov) paste0(
"  for (q in 1:Q) {
    if (x_binary[q] == 1) {
      beta_f[q] ~ ", pr$beta_binary, ";
    } else {
      beta_f[q] ~ ", pr$beta, ";
    }
  }
") else "",
if (has_cov && has_hyper) paste0(
"  for (q in 1:Q) {
    if (x_binary[q] == 1) {
      beta_g[q] ~ ", pr$beta_binary, ";
    } else {
      beta_g[q] ~ ", pr$beta, ";
    }
  }
") else "",
if (identical(likelihood, "ordinal")) paste0(
"  for (p in 1:P) cut[p] ~ ", pr$cut, ";
") else "",
likelihood_block,
"}
")

  gq_block <- paste0(
"generated quantities {
  // In-iterations sign relabeling (blavaan convention; Merkle et al.,
  // 2021; evaluated in Chen, Miocevic, & Falk, 2025): the likelihood is
  // reflection-invariant, so raw *_un draws may occupy mirror modes
  // across chains BY DESIGN. Per iteration, each factor is reflected so
  // its reference item loading is positive, with all sign-bearing
  // quantities of that factor reflected jointly; hyper orientations
  // follow their reference child. Report and diagnose these corrected
  // quantities, never the *_un scaffolding.
  vector[K] sign_f;
  vector[P] lambda;
",
if (has_cov)
"  matrix[Q, K] beta_f_cor;
" else "",
if (has_hyper)
"  vector[G] sign_g;
" else "",
if (free_hyper)
"  vector[K] delta;
  vector[K] delta_std;
" else "",
if (has_cov && has_hyper)
"  matrix[Q, G] beta_g_cor;
" else "",
if (has_hyper && !free_hyper)
"  vector[G] rho_g;
" else "",
"  for (k in 1:K) {
    sign_f[k] = lambda_un[ref_item[k]] < 0 ? -1 : 1;
  }
  for (p in 1:P) {
    lambda[p] = sign_f[item_factor[p]] * lambda_un[p];
  }
",
if (has_cov)
"  for (k in 1:K) {
    beta_f_cor[, k] = sign_f[k] * beta_f[, k];
  }
" else "",
if (free_hyper)
"  {
    vector[K] delta_tmp;
    for (k in 1:K) delta_tmp[k] = sign_f[k] * delta_l[k];
    for (g in 1:G) sign_g[g] = delta_tmp[ref_child[g]] < 0 ? -1 : 1;
    for (k in 1:K) {
      delta[k] = factor_hyper[k] == 0 ? 0
                 : sign_g[factor_hyper[k]] * delta_tmp[k];
      delta_std[k] = delta[k] / sqrt(square(delta[k]) + 1);
    }
  }
" else "",
if (has_hyper && !free_hyper)
"  // Unity loadings: per-factor reflection is only a joint block
  // symmetry, so the hyper orientation follows its reference child.
  for (g in 1:G) sign_g[g] = sign_f[ref_child[g]];
  for (g in 1:G) rho_g[g] = square(tau_g[g]) / (square(tau_g[g]) + 1);
" else "",
if (has_cov && has_hyper)
"  for (g in 1:G) {
    beta_g_cor[, g] = sign_g[g] * beta_g[, g];
  }
" else "",
"}
")

  paste0(data_block, params_block, transformed_block, model_block, gq_block)
}

# Augment the base Stan data list for the GPU-vectorized program:
# precompute the linear gather index into to_vector(eta) (column-major),
# and for the ordinal likelihood sort the long format by item and record
# contiguous per-item segments (so each ordered_logistic call shares one
# cutpoint vector). This is the automatic vectorization/matrix-layout
# transformation applied whenever opencl_ids is supplied.
.dimsem_stan_data_gpu <- function(sd, likelihood = "normal") {
  # Column-major flattening: to_vector(eta)[(k - 1) * N + r] == eta[r, k].
  fac_of_obs <- sd$item_factor[sd$obs_col]
  sd$eta_idx <- (fac_of_obs - 1L) * sd$N + sd$obs_row

  if (identical(likelihood, "ordinal")) {
    ord <- order(sd$obs_col)
    sd$obs_row <- sd$obs_row[ord]
    sd$obs_col <- sd$obs_col[ord]
    sd$eta_idx <- sd$eta_idx[ord]
    sd$y_int <- as.integer(sd$y_obs[ord])
    sd$y_obs <- NULL
    seg_start <- integer(sd$P)
    seg_end <- integer(sd$P)
    for (p in seq_len(sd$P)) {
      w <- which(sd$obs_col == p)
      if (length(w) > 0) {
        seg_start[p] <- w[1]
        seg_end[p] <- w[length(w)]
      } else {
        seg_start[p] <- 1L
        seg_end[p] <- 0L
      }
    }
    sd$seg_start <- seg_start
    sd$seg_end <- seg_end
  }

  sd$grainsize <- NULL
  sd
}


# Delete the cached executable(s) belonging to a Stan source file so the
# next cmdstan_model() call performs a full rebuild. cmdstanr places the
# binary next to the (hash-named) .stan file: same path without the
# extension on Unix-alikes, plus an .exe suffix on Windows. Pure file
# logic by design -- no cmdstanr dependency -- so it can be unit-tested
# and reused.
.dimsem_clear_stan_exe <- function(stan_file, verbose = TRUE) {
  exe <- sub("\\.stan$", "", stan_file)
  candidates <- unique(c(exe, paste0(exe, ".exe")))
  removed <- candidates[file.exists(candidates)]
  if (length(removed) > 0) {
    unlink(removed)
    if (isTRUE(verbose)) {
      message("force_recompile = TRUE: removed cached Stan executable (",
              paste(basename(removed), collapse = ", "),
              "); rebuilding with the current compile options.")
    }
  } else if (isTRUE(verbose)) {
    message("force_recompile = TRUE: no cached Stan executable found; ",
            "compiling fresh.")
  }
  invisible(removed)
}


# --- Posterior summaries for Bayesian DimSem fits ---------------------------------

# Build a lookup translating Stan parameter indices back to SEM labels, so
# posterior summaries read in item/factor/covariate terms rather than raw
# subscripts. Depends only on the estimation spec and covariate set, so it
# is cheap to store on the fitted object.
.dimsem_stan_param_map <- function(spec, covariates = NULL,
                                   hyper_loadings = "unity") {
  items <- names(spec$partition)[!is.na(spec$partition)]
  part <- spec$partition[items]
  latent_labels <- paste0("F", sort(unique(part)))
  hyper_names <- names(spec$hyper_blocks)

  list(
    items = items,
    factors = latent_labels,
    hyper = hyper_names,
    covariates = covariates %||% character(0),
    hyper_loadings = hyper_loadings
  )
}

#' Extract core posterior summaries from a Bayesian DimSem fit
#'
#' Relabels the posterior summary of a [DimSEM_estimate()] Bayesian fit into
#' interpretable SEM terms (loadings as `F =~ item`, structural paths as
#' `factor ~ covariate`, hyper loadings/correlations), each with posterior
#' mean/median/sd, a central credible interval, and convergence
#' diagnostics (rhat, ess_bulk, ess_tail).
#'
#' @param object A `"DimSem_estimate"` object fitted with `engine = "bayes"`.
#' @param pars `"core"` (default: loadings, structural, hyper) or `"all"`
#'   (adds item intercepts and residual SDs).
#' @param prob Central credible-interval mass (default .95).
#' @return A `"DimSEM_posterior"` list of tidy data frames, one per block.
#' @export
# Extract core posterior summary statistics from a DimSem_estimate object,
# relabeled into SEM terms. Returns a list of tidy data frames by block
# (loadings, intercepts, residual_sd, structural, hyper, correlations),
# each carrying posterior mean/median/sd, a central credible interval, and
# convergence diagnostics (rhat, ess_bulk, ess_tail).
#
#   object : a "DimSem_estimate" object fitted with engine = "bayes".
#   pars   : which blocks to return; "core" (default) = loadings,
#            structural, hyper, correlations; "all" adds intercepts and
#            residual SDs; or a character vector of block names.
#   prob   : central credible-interval mass (default .95).
DimSEM_posterior <- function(object, pars = c("core", "all"), prob = 0.95) {
  if (!inherits(object, "DimSem_estimate") &&
      !inherits(object, "DimSEM_estimate")) {
    stop("`object` must be a DimSem_estimate object.", call. = FALSE)
  }
  if (!identical(object$engine, "bayes")) {
    stop("Posterior summaries are only available for engine = \"bayes\"; ",
         "this fit used engine = \"", object$engine, "\". Use ",
         "summary(object$fit) / lavaan tools for ML fits.", call. = FALSE)
  }
  fit <- object$fit
  if (is.null(fit)) {
    stop("No fitted Stan object is present (a dry run stores only the ",
         "generated program in `$stan`). Re-run with a real cmdstanr ",
         "backend to obtain posterior draws.", call. = FALSE)
  }
  pars <- if (length(pars) > 1 && all(pars == c("core", "all"))) {
    "core"
  } else {
    match.arg(pars[1], c("core", "all"))
  }
  map <- object$param_map
  if (is.null(map)) {
    stop("The fitted object lacks a parameter map; re-fit with a DimSem ",
         "version that stores `param_map`.", call. = FALSE)
  }

  tail_p <- (1 - prob) / 2
  qs <- c(tail_p, 1 - tail_p)
  # cmdstanr's summary() takes posterior-style quantile functions.
  summ <- function(vars) {
    s <- tryCatch(
      fit$summary(variables = vars, "mean", "median", "sd",
                  q_lo = ~ stats::quantile(.x, probs = qs[1]),
                  q_hi = ~ stats::quantile(.x, probs = qs[2]),
                  "rhat", "ess_bulk", "ess_tail"),
      error = function(e) NULL
    )
    if (is.null(s)) NULL else as.data.frame(s)
  }

  ci_names <- c(paste0("q", round(100 * qs[1], 1)),
                paste0("q", round(100 * qs[2], 1)))
  relabel <- function(df, labels) {
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df$variable <- labels
    names(df)[names(df) %in% c("q_lo", "q_hi")] <- ci_names
    df
  }

  out <- list()

  # Loadings: lambda[p] -> "F<k> =~ <item>".
  part <- object$proposal_partition
  lam <- summ("lambda")
  if (!is.null(lam)) {
    lab <- vapply(seq_along(map$items), function(i) {
      paste0("F", part[map$items[i]], " =~ ", map$items[i])
    }, character(1))
    out$loadings <- relabel(lam, lab)
  }

  if (identical(pars, "all")) {
    nu <- summ("nu")
    if (!is.null(nu)) out$intercepts <- relabel(nu, paste0(map$items, " ~1"))
    sig <- summ("sigma")
    if (!is.null(sig)) {
      out$residual_sd <- relabel(sig, paste0(map$items, " ~~ ", map$items))
    }
  }

  # Structural coefficients: beta_f[q, k] -> "<factor> ~ <covariate>".
  if (length(map$covariates) > 0) {
    bf <- summ("beta_f_cor")
    if (!is.null(bf)) {
      grid <- expand.grid(covariate = map$covariates, factor = map$factors,
                          stringsAsFactors = FALSE)
      out$structural <- relabel(bf, paste0(grid$factor, " ~ ", grid$covariate))
    }
    if (length(map$hyper) > 0) {
      bg <- summ("beta_g_cor")
      if (!is.null(bg)) {
        gridg <- expand.grid(covariate = map$covariates, hyper = map$hyper,
                             stringsAsFactors = FALSE)
        out$structural_hyper <- relabel(
          bg, paste0(gridg$hyper, " ~ ", gridg$covariate))
      }
    }
  }

  # Hyper level: unity -> rho_g (implied correlation); free -> delta_std.
  if (length(map$hyper) > 0) {
    if (identical(map$hyper_loadings, "unity")) {
      rg <- summ("rho_g")
      if (!is.null(rg)) {
        out$hyper_correlation <- relabel(
          rg, paste0(map$hyper, " implied within-block r"))
      }
    } else {
      ds <- summ("delta_std")
      if (!is.null(ds)) {
        # One standardized hyper loading per first-order factor.
        blocks <- object$hyper_blocks %||% list()
        lab <- vapply(map$factors, function(fk) {
          owner <- names(which(vapply(blocks, function(b) fk %in% b,
                                      logical(1))))
          if (length(owner) == 0) paste0(fk, " (stand-alone)")
          else paste0(owner[1], " =~ ", fk)
        }, character(1))
        out$hyper_loadings_std <- relabel(ds, lab)
      }
    }
  }

  attr(out, "prob") <- prob
  attr(out, "engine") <- object$engine
  class(out) <- "DimSEM_posterior"
  out
}

#' @export
print.DimSEM_posterior <- function(x, ...) {
  cat("DimSem posterior summaries (",
      round(100 * attr(x, "prob")), "% credible intervals)\n\n", sep = "")
  for (block in names(x)) {
    if (is.null(x[[block]]) || nrow(x[[block]]) == 0) next
    cat("## ", block, "\n", sep = "")
    df <- x[[block]]
    num <- vapply(df, is.numeric, logical(1))
    df[num] <- lapply(df[num], function(v) round(v, 3))
    print(df, row.names = FALSE)
    cat("\n")
  }
  invisible(x)
}

# Choose the sign-reference item for each factor: by default the marginal
# hub (highest mean |Pearson r| with its own community members), i.e. the
# item most reliably far from a zero loading; overridable via a named
# vector sign_ref = c(F1 = "item_a", G1 = "F2", ...).
.dimsem_sign_refs <- function(data, items, factor_index, sign_ref = NULL) {
  K <- max(factor_index)
  R <- suppressWarnings(stats::cor(as.matrix(data[, items, drop = FALSE]),
                                   use = "pairwise.complete.obs"))
  R[!is.finite(R)] <- 0
  ref <- integer(K)
  for (k in seq_len(K)) {
    members <- which(factor_index == k)
    if (length(members) == 1L) {
      ref[k] <- members
      next
    }
    hub <- vapply(members, function(i) {
      mean(abs(R[items[i], items[setdiff(members, i)]]))
    }, numeric(1))
    ref[k] <- members[which.max(hub)]
    lab <- paste0("F", k)
    if (!is.null(sign_ref) && lab %in% names(sign_ref)) {
      cand <- match(sign_ref[[lab]], items)
      if (!is.na(cand) && cand %in% members) ref[k] <- cand
    }
  }
  ref
}

#' Lavaan-style parameter table for Bayesian DimSem fits
#'
#' Returns the main model parameters of a [DimSEM_estimate()] Bayesian fit
#' in the layout of [lavaan::parameterEstimates()]: `lhs`/`op`/`rhs` rows
#' (`=~` loadings including hyper loadings, `~1` intercepts, `~~` residual
#' and hyper-factor variances, `~` covariate regressions, `:=` defined
#' quantities such as the implied within-block correlation under unity
#' loadings), excluding all Stan scaffolding (unconstrained `*_un`
#' parameters, latent scores, innovations, sign indicators). Sign-bearing
#' quantities come from the sign-relabeled generated quantities, so signs
#' and diagnostics are directly interpretable. Columns: posterior mean
#' (`est`), posterior SD (`se`), central credible interval at `prob`, and
#' `rhat`/`ess_bulk`/`ess_tail`. Fixed values (unity hyper loadings,
#' standardized disturbances) show `se = 0` and `NA` diagnostics,
#' mirroring lavaan's display of fixed parameters.
#'
#' @param object A `"DimSem_estimate"` object fitted with `engine = "bayes"`.
#' @param prob Central credible-interval mass (default .95).
#' @return A data frame in lavaan `parameterEstimates()` layout.
#' @export
DimSEM_parameterEstimates <- function(object, prob = 0.95) {
  if (!inherits(object, "DimSem_estimate") &&
      !inherits(object, "DimSEM_estimate")) {
    stop("`object` must be a DimSem_estimate object.", call. = FALSE)
  }
  if (!identical(object$engine, "bayes")) {
    stop("Use lavaan::parameterEstimates(object$fit) for ML fits; this ",
         "table is for engine = \"bayes\".", call. = FALSE)
  }
  if (is.null(object$fit)) {
    stop("No fitted Stan object present (dry run); re-run with a real ",
         "cmdstanr backend.", call. = FALSE)
  }
  .dimsem_require_namespace("posterior")
  map <- object$param_map
  part <- object$proposal_partition
  blocks <- object$hyper_blocks %||% list()
  unity <- identical(map$hyper_loadings, "unity")
  qs <- c((1 - prob) / 2, 1 - (1 - prob) / 2)

  summarize_var <- function(a) {
    m <- as.numeric(a)
    c(est = mean(m), se = stats::sd(m),
      ci.lower = unname(stats::quantile(m, qs[1])),
      ci.upper = unname(stats::quantile(m, qs[2])),
      rhat = posterior::rhat(posterior::as_draws_array(a)),
      ess_bulk = posterior::ess_bulk(posterior::as_draws_array(a)),
      ess_tail = posterior::ess_tail(posterior::as_draws_array(a)))
  }
  grab <- function(vars) {
    d <- tryCatch(object$fit$draws(variables = vars),
                  error = function(e) NULL)
    if (is.null(d)) NULL else posterior::as_draws_array(d)
  }
  rows <- list()
  add <- function(lhs, op, rhs, stats_row) {
    rows[[length(rows) + 1]] <<- data.frame(
      lhs = lhs, op = op, rhs = rhs, t(stats_row),
      stringsAsFactors = FALSE)
  }
  add_fixed <- function(lhs, op, rhs, value) {
    add(lhs, op, rhs, c(est = value, se = 0, ci.lower = value,
                        ci.upper = value, rhat = NA_real_,
                        ess_bulk = NA_real_, ess_tail = NA_real_))
  }

  lam <- grab("lambda")
  if (!is.null(lam)) {
    for (i in seq_along(map$items)) {
      add(paste0("F", part[map$items[i]]), "=~", map$items[i],
          summarize_var(lam[, , i, drop = FALSE]))
    }
  }
  if (length(blocks) > 0) {
    if (unity) {
      for (g in names(blocks)) {
        for (fk in blocks[[g]]) add_fixed(g, "=~", fk, 1)
      }
    } else {
      ds <- grab("delta")
      if (!is.null(ds)) {
        for (g in names(blocks)) {
          for (fk in blocks[[g]]) {
            k <- match(fk, map$factors)
            add(g, "=~", fk, summarize_var(ds[, , k, drop = FALSE]))
          }
        }
      }
    }
  }
  nu <- grab("nu")
  if (!is.null(nu)) {
    for (i in seq_along(map$items)) {
      add(map$items[i], "~1", "", summarize_var(nu[, , i, drop = FALSE]))
    }
  }
  # Residual VARIANCES: sigma^2 transformed draw-wise (lavaan-canonical).
  sig <- grab("sigma")
  if (!is.null(sig)) {
    for (i in seq_along(map$items)) {
      add(map$items[i], "~~", map$items[i],
          summarize_var(sig[, , i, drop = FALSE]^2))
    }
  }
  for (fk in map$factors) add_fixed(fk, "~~", fk, 1)
  if (length(blocks) > 0) {
    if (unity) {
      tg <- grab("tau_g")
      if (!is.null(tg)) {
        for (gi in seq_along(names(blocks))) {
          add(names(blocks)[gi], "~~", names(blocks)[gi],
              summarize_var(tg[, , gi, drop = FALSE]^2))
        }
      }
    } else {
      for (g in names(blocks)) add_fixed(g, "~~", g, 1)
    }
  }
  if (length(map$covariates) > 0) {
    bf <- grab("beta_f_cor")
    if (!is.null(bf)) {
      for (k in seq_along(map$factors)) {
        for (q in seq_along(map$covariates)) {
          idx <- match(paste0("beta_f_cor[", q, ",", k, "]"),
                       dimnames(bf)[[3]])
          if (!is.na(idx)) {
            add(map$factors[k], "~", map$covariates[q],
                summarize_var(bf[, , idx, drop = FALSE]))
          }
        }
      }
    }
    if (length(blocks) > 0) {
      bg <- grab("beta_g_cor")
      if (!is.null(bg)) {
        for (gi in seq_along(names(blocks))) {
          for (q in seq_along(map$covariates)) {
            idx <- match(paste0("beta_g_cor[", q, ",", gi, "]"),
                         dimnames(bg)[[3]])
            if (!is.na(idx)) {
              add(names(blocks)[gi], "~", map$covariates[q],
                  summarize_var(bg[, , idx, drop = FALSE]))
            }
          }
        }
      }
    }
  }
  if (length(blocks) > 0 && unity) {
    rg <- grab("rho_g")
    if (!is.null(rg)) {
      for (gi in seq_along(names(blocks))) {
        add(paste0("rho_", names(blocks)[gi]), ":=", "",
            summarize_var(rg[, , gi, drop = FALSE]))
      }
    }
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  attr(out, "prob") <- prob
  out
}
