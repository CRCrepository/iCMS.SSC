#' Classify a single sample using the K-Nearest Neighbours (KNN) method
#'
#' @param v Named numeric vector of (log-transformed) gene expression
#'     values.
#' @param nn Number of nearest neighbours (default 10).
#' @param metric Similarity metric: "pearson", "kendall" (default),
#'     "spearman", or "cosine".
#' @param allgenes Use extended gene set (default FALSE; uses
#'     epithelial genes only).
#' @param min.genes Minimum gene overlap with the prototype set
#'     (default 30).
#' @param mcenter Mean-centre the sample before computing
#'     similarity. Only relevant for cosine similarity (correlation
#'     metrics are shift-invariant). Default FALSE.
#' @param stratify.ds If TRUE, vote is stratified by training dataset:
#'     for each dataset the class whose best prototype is most similar
#'     wins one vote, then votes are tallied across datasets. This
#'     prevents a single batch-matched dataset from dominating the
#'     top-nn neighbours. You may consider this if the input data have
#'     strong batch effects. Default FALSE.
#' @param weighted.genes If TRUE, restrict computation to the top 100
#'     most discriminative genes (by absolute t-statistic); for
#'     "pearson" and "cosine" uses weighted correlation / weighted
#'     cosine instead. Not useful if you are using the recommended
#'     "kendall" distance metric. Default FALSE.
#' @param verbose If TRUE, return individual neighbour labels,
#'     similarities, and dataset labels as additional columns. Default
#'     FALSE.
#' @return A one-row data frame with columns: qi2, qi3, i2, i3, i2i3,
#'     nearest.icms, confident.icms, nearest.ds, ngenes (plus
#'     neighbour details when verbose = TRUE).
#' @export
knn.vect <- function(v, nn, metric, allgenes, min.genes,
                     mcenter = FALSE, stratify.ds = FALSE,
                     weighted.genes = FALSE, verbose = FALSE) {

  prot <- if (allgenes) ssproto.ext else ssproto
  cg   <- intersect(rownames(prot), names(na.omit(v)))

  if (length(cg) < min.genes) {
    warning("Too few common genes. Unable to map prototype.")
    return(NA)
  }

  if (mcenter) v <- v - mean(v, na.rm = TRUE)

  ## --- Gene weighting / subsetting -----------------------------------
  if (weighted.genes) {
    wts <- gene.weights[cg]
    if (metric %in% c("pearson", "cosine")) {
      smat_vec <- if (metric == "pearson")
        .wt_pearson_mat(v[cg], prot[cg, ], wts)
      else
        .wt_cosine_mat(v[cg], prot[cg, ], wts)
      smat <- matrix(smat_vec, nrow = 1, dimnames = list(NULL, colnames(prot)))
    } else {
      top_cg <- head(cg[order(wts, decreasing = TRUE)], min(100L, length(cg)))
      if (length(top_cg) < min.genes) {
        warning("Too few genes after weighting. Unable to map prototype.")
        return(NA)
      }
      cg   <- top_cg
      smat <- cor(as.matrix(v)[cg, ], prot[cg, ], method = metric)
      colnames(smat) <- colnames(prot)
    }
  } else {
    if (metric != "cosine") {
      smat <- cor(as.matrix(v)[cg, ], prot[cg, ], method = metric)
      colnames(smat) <- colnames(prot)
    } else {
      smat <- .cosine_vec_mat(as.matrix(v)[cg, ], prot[cg, ])
      colnames(smat) <- colnames(prot)
    }
  }

  ## Dataset labels for each prototype column
  ds_all <- sub("^i[23][.]", "", sub("-s[0-9]+$", "", colnames(prot)))

  ## Stratify voting by dataset
  if (stratify.ds) {
    unique_ds  <- unique(ds_all)
    ds_votes   <- vapply(unique_ds, function(d) {
      mask  <- ds_all == d
      sims  <- as.vector(smat)[mask]
      idx_d <- ssicms.index[mask]
      if (max(sims[idx_d == "i2"]) > max(sims[idx_d == "i3"])) "i2" else "i3"
    }, character(1L))
    names(ds_votes) <- unique_ds

    n_ds  <- length(unique_ds)
    n_i2  <- sum(ds_votes == "i2")
    n_i3  <- sum(ds_votes == "i3")

    ## Weighted proportions (using per-dataset max similarity as weights)
    ds_strength <- vapply(unique_ds, function(d) {
      mask  <- ds_all == d
      sims  <- as.vector(smat)[mask]
      max(sims)
    }, numeric(1L))
    wt_i2 <- sum(ds_strength[ds_votes == "i2"]) / sum(ds_strength)
    wt_i3 <- sum(ds_strength[ds_votes == "i3"]) / sum(ds_strength)

    nearest.icms   <- if (n_i2 >= n_i3) "i2" else "i3"
    ## Confident when unanimous (all training datasets agree)
    confident.icms <- if (n_i2 == n_ds || n_i3 == n_ds) nearest.icms else NA

    ## Summary similarities: median of per-dataset max similarities
    qi2 <- median(vapply(unique_ds, function(d) {
      max(as.vector(smat)[ds_all == d & ssicms.index == "i2"])
    }, numeric(1L)))
    qi3 <- median(vapply(unique_ds, function(d) {
      max(as.vector(smat)[ds_all == d & ssicms.index == "i3"])
    }, numeric(1L)))

    nearest.ds <- unique_ds[which.max(ds_strength)]
    i2i3       <- wt_i2 - wt_i3

    retl <- list(
      qi2            = qi2,
      qi3            = qi3,
      i2             = wt_i2,
      i3             = wt_i3,
      i2i3           = i2i3,
      nearest.icms   = nearest.icms,
      confident.icms = confident.icms,
      nearest.ds     = nearest.ds,
      ngenes         = length(cg)
    )
    return(data.frame(retl, stringsAsFactors = FALSE))
  }

  ## Without dataset stratification (default)
  ord     <- order(smat, decreasing = TRUE)[1:nn]
  knn_cls <- ssicms.index[ord]
  knnd    <- as.vector(smat[ord])
  names(knnd) <- colnames(smat)[ord]
  knndset <- ds_all[ord]

  ## Distance-weighted vote: weight each neighbour by its similarity score.
  ## Shift similarities to be non-negative so weights are valid.
  w_shift <- pmax(knnd - min(0, min(knnd)), 0)
  w_total <- sum(w_shift)
  if (w_total > 0) {
    wt_i2 <- sum(w_shift[knn_cls == "i2"]) / w_total
    wt_i3 <- sum(w_shift[knn_cls == "i3"]) / w_total
  } else {
    wt_i2 <- wt_i3 <- 0.5
  }

  ## Raw counts (used for binomial confidence test)
  n_i2 <- sum(knn_cls == "i2")
  n_i3 <- sum(knn_cls == "i3")

  ## Binomial confidence: consistent across different nn values.
  ## H0: p(i2) = 0.5. Confident when one-tailed p < 0.05.
  n_majority <- max(n_i2, n_i3)
  p_binom    <- pbinom(n_majority - 1L, nn, 0.5, lower.tail = FALSE)
  confident.icms <- if (p_binom < 0.05)
    if (n_i2 > n_i3) "i2" else "i3"
  else NA

  nearest.icms <- if (wt_i2 >= wt_i3) "i2" else "i3"
  i2i3         <- wt_i2 - wt_i3

  ## Median similarity per class (summary of neighbour quality)
  qd  <- tapply(knnd, knn_cls, median)
  qi2 <- if ("i2" %in% names(qd)) as.vector(qd["i2"]) else NA_real_
  qi3 <- if ("i3" %in% names(qd)) as.vector(qd["i3"]) else NA_real_

  if (verbose) {
    retl <- list(nn = knn_cls, nnd = knnd, nnds = knndset)
  } else {
    retl <- list()
  }

  retl <- c(retl,
    qi2            = qi2,
    qi3            = qi3,
    i2             = wt_i2,
    i3             = wt_i3,
    i2i3           = i2i3,
    nearest.icms   = nearest.icms,
    confident.icms = confident.icms,
    nearest.ds     = knndset[1],
    ngenes         = length(cg)
  )
  data.frame(retl, stringsAsFactors = FALSE)
}
