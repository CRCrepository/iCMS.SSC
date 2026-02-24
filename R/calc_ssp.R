## Internal helper: cosine similarity between a single vector x and every column of matrix y.
## x : numeric vector (length p)
## y : numeric matrix (p x n)
## Returns a 1 x n matrix.
.cosine_vec_mat <- function(x, y) {
  a  <- crossprod(x, y)          # 1 x n dot products
  nx <- norm(x, type = "2")
  ny <- apply(y, 2, norm, type = "2")
  a / (nx * ny)
}

## Internal helper: weighted Pearson correlation between vector x and every column of matrix y.
## w : positive numeric weights (length p), will be normalised internally.
.wt_pearson_mat <- function(x, y, w) {
  w   <- w / sum(w)
  xm  <- sum(w * x)
  ym  <- colSums(w * y)              # weighted column means
  xc  <- x - xm
  yc  <- sweep(y, 2, ym, "-")        # centre each prototype column
  cov <- as.vector(crossprod(w * xc, yc))
  vx  <- sum(w * xc^2)
  vy  <- colSums(w * yc^2)
  cov / sqrt(vx * vy)
}

## Internal helper: weighted cosine similarity between vector x and every column of matrix y.
.wt_cosine_mat <- function(x, y, w) {
  ws <- sqrt(w / sum(w))
  xw <- x * ws
  yw <- y * ws                       # broadcast weights to each row
  nv <- sqrt(sum(xw^2))
  ny <- sqrt(colSums(yw^2))
  as.vector(crossprod(xw, yw)) / (nv * ny)
}


#' Classify a single sample using the Distance-Quantile (DQ) single-sample predictor
#'
#' @param v Named numeric vector of (log-transformed) gene expression values.
#' @param min.cor Minimum q-quantile similarity required for a confident call (default 0.1).
#' @param min.dist Minimum scaled margin between i2 and i3 quantile scores for a confident call (default 0.05).
#' @param q Quantile used to summarise per-class similarity scores (default 0.9).
#'   Values 0.9–1 suit non-batch-corrected data; 0.6–0.7 suit batch-corrected data.
#' @param metric Similarity metric: "pearson", "kendall" (default), "spearman", or "cosine".
#' @param allgenes Use extended gene set (default FALSE; uses epithelial genes only).
#' @param min.genes Minimum gene overlap with the prototype set (default 30).
#' @param stratify.ds If TRUE, compute the quantile margin per training dataset and aggregate
#'   across datasets, which removes the confounding effect of platform-level similarity
#'   differences (default FALSE).
#' @param weighted.genes If TRUE, restrict computation to the top 100 most discriminative genes
#'   (by absolute t-statistic); for "pearson" and "cosine" uses weighted correlation /
#'   weighted cosine instead (default FALSE).
#' @param verbose If TRUE, return all 192 prototype-level similarity columns in addition to
#'   the summary statistics. Default FALSE returns only the 7 summary columns.
#' @return A one-row data frame with columns: i2, i3, i2i3, nearest.icms, confident.icms,
#'   nearest.dset, ngenes (plus prototype columns when verbose = TRUE).
#' @export
dq.vect <- function(v, min.cor, min.dist, q, metric, allgenes, min.genes,
                    stratify.ds = FALSE, weighted.genes = FALSE, verbose = FALSE) {

  prot <- if (allgenes) ssproto.ext else ssproto
  cg   <- intersect(rownames(prot), names(na.omit(v)))

  if (length(cg) < min.genes) {
    warning("Too few common genes. Unable to map prototype.")
    return(NA)
  }

  ## --- Gene weighting / subsetting -----------------------------------
  if (weighted.genes) {
    wts <- gene.weights[cg]               # named vector of |t-statistic|
    if (metric %in% c("pearson", "cosine")) {
      ## Use weighted similarity (all genes, but discriminative ones count more)
      if (metric == "pearson") {
        smat_vec <- .wt_pearson_mat(v[cg], prot[cg, ], wts)
      } else {
        smat_vec <- .wt_cosine_mat(v[cg], prot[cg, ], wts)
      }
      smat <- matrix(smat_vec, nrow = 1, dimnames = list(NULL, colnames(prot)))
    } else {
      ## For rank-based metrics: restrict to top 100 most discriminative genes
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
    ## Standard unweighted similarity
    if (metric != "cosine") {
      smat <- cor(as.matrix(v)[cg, ], prot[cg, ], method = metric)
      colnames(smat) <- colnames(prot)
    } else {
      smat <- .cosine_vec_mat(as.matrix(v)[cg, ], prot[cg, ])
      colnames(smat) <- colnames(prot)
    }
  }

  ## --- Dataset labels ------------------------------------------------
  ds_all <- sub("^i[23][.]", "", sub("-s[0-9]+$", "", colnames(prot)))

  ## --- Per-class quantile summaries ----------------------------------
  if (stratify.ds) {
    ## Compute quantile margin separately for each training dataset,
    ## then aggregate — removes dataset-level similarity lift.
    unique_ds <- unique(ds_all)
    margins <- vapply(unique_ds, function(d) {
      sims_d <- as.vector(smat)[ds_all == d]
      idx_d  <- ssicms.index[ds_all == d]
      qi2_d  <- quantile(sims_d[idx_d == "i2"], q, na.rm = TRUE)
      qi3_d  <- quantile(sims_d[idx_d == "i3"], q, na.rm = TRUE)
      qi2_d - qi3_d
    }, numeric(1L))
    names(margins) <- unique_ds

    score <- mean(margins)
    qi2   <- mean(vapply(unique_ds, function(d) {
      sims_d <- as.vector(smat)[ds_all == d & ssicms.index == "i2"]
      quantile(sims_d, q, na.rm = TRUE)
    }, numeric(1L)))
    qi3   <- mean(vapply(unique_ds, function(d) {
      sims_d <- as.vector(smat)[ds_all == d & ssicms.index == "i3"]
      quantile(sims_d, q, na.rm = TRUE)
    }, numeric(1L)))
    i2i3  <- score

    nearest.icms  <- if (score >= 0) "i2" else "i3"
    ## Confident when all datasets agree AND mean absolute signal meets threshold
    all_agree <- all(margins > 0) || all(margins < 0)
    confident.icms <- if (all_agree && max(qi2, qi3) >= min.cor &&
                          sqrt(2 * score^2) >= min.dist) nearest.icms else NA

    nearest.dset <- unique_ds[which.max(abs(margins))]  # dataset with strongest signal

  } else {
    ## Global quantile (original DQ logic)
    maxdist <- t(apply(smat, 1, function(x)
      tapply(x, ssicms.index, quantile, q, na.rm = TRUE)))
    qi2   <- maxdist[1, "i2"]
    qi3   <- maxdist[1, "i3"]
    i2i3  <- qi2 - qi3

    nearest.icms <- if (qi2 >= qi3) "i2" else "i3"

    ## Confidence criteria (applied sequentially):
    ## (1) Negative correlation with both classes → no call (sample has no expression signal)
    ## (2) Both: very weak similarity AND very small margin → no call
    ##     Note: a sample with max < min.cor but a large margin still gets a call;
    ##     this correctly handles sparse data where absolute correlations are low but
    ##     directional signal is present.
    ## (3) Scaled margin too small → no call
    ##     dd = sqrt(2) * |i2 - i3|, the distance from the i2=i3 diagonal scaled by sqrt(2)
    dd <- sqrt(2 * (qi2 - qi3)^2)
    confident.icms <- if (max(qi2, qi3) < 0 ||
                          (max(qi2, qi3) < min.cor & abs(qi2 - qi3) < min.cor) ||
                          dd < min.dist) NA else nearest.icms

    nearest.dset <- ds_all[which.max(as.vector(smat))]
  }

  ## --- Build output ---------------------------------------------------
  summary_cols <- data.frame(
    i2             = qi2,
    i3             = qi3,
    i2i3           = i2i3,
    nearest.icms   = nearest.icms,
    confident.icms = confident.icms,
    nearest.dset   = nearest.dset,
    ngenes         = length(cg),
    stringsAsFactors = FALSE
  )

  if (verbose) {
    proto_cols <- as.data.frame(smat)
    return(cbind(proto_cols, summary_cols))
  }
  summary_cols
}


#' Calculate cosine similarity between two matrices
#'
#' @param x Numeric matrix (genes x samples).
#' @param y Numeric matrix (genes x references).
#' @return Cosine similarity matrix (samples x references).
#' @export
cosine.sim <- function(x, y) {
  a <- crossprod(x, y)
  b <- outer(sqrt(colSums(x^2)), sqrt(colSums(y^2)))
  a / b
}
