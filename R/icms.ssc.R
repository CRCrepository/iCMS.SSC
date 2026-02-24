#' Classify bulk gene-expression samples using the Distance-Quantile (DQ) single-sample predictor
#'
#' Estimates the most probable iCMS subtype for each sample by comparing its similarity to
#' synthetic exemplar centroids from three training cohorts (Marisa, PETACC, TCGA), using
#' the q-quantile of per-class similarities as the decision score.
#'
#' @param ivect Gene-expression matrix with gene symbols as row names and samples as columns.
#'   Values should be log-transformed.
#' @param min.cor Minimum quantile similarity required for a confident call (default 0.1).
#' @param min.dist Minimum scaled margin between i2 and i3 scores for a confident call (default 0.05).
#' @param q Quantile used to summarise per-class similarities (default 0.9).
#'   Values 0.9–1 suit non-batch-corrected data; 0.6–0.7 suit batch-corrected data.
#' @param metric Similarity metric: "pearson", "kendall" (default), "spearman", or "cosine".
#' @param allgenes Use extended gene set (default FALSE; uses epithelial genes only).
#' @param min.genes Minimum gene overlap required (default 30).
#' @param stratify.ds If TRUE, compute per-dataset quantile margins and aggregate, removing
#'   the confounding effect of platform-level similarity differences (default FALSE).
#' @param weighted.genes If TRUE, up-weight discriminative genes by their absolute
#'   t-statistic when metric is "pearson" or "cosine"; for rank-based metrics restricts
#'   to the top 100 most discriminative genes (default FALSE).
#' @param jobs Number of parallel cores. Set to 1 to disable parallelism (default 4).
#' @param verbose If TRUE, include all 192 prototype-level similarity columns in the output.
#'   Default FALSE returns 7 summary columns only.
#' @return A data frame with one row per sample. Key columns: nearest.icms, confident.icms.
#' @export
iCMS.DQ <- function(ivect,
                    min.cor       = 0.1,
                    min.dist      = 0.05,
                    q             = 0.9,
                    metric        = "kendall",
                    allgenes      = FALSE,
                    min.genes     = 30,
                    stratify.ds   = FALSE,
                    weighted.genes = FALSE,
                    jobs          = 4,
                    verbose       = FALSE) {

  message("input matrix: ", ncol(ivect), " samples, ", nrow(ivect),
          " genes. DQ method, metric=", metric,
          if (stratify.ds) ", stratify.ds=TRUE" else "",
          if (weighted.genes) ", weighted.genes=TRUE" else "")

  if (q > 1) q <- 1
  if (q < 0) q <- 0
  if (q < 0.5) warning("q < 0.5: considering similarities below the median.")

  valid.metrics <- c("pearson", "kendall", "spearman", "cosine")
  if (!metric %in% valid.metrics) {
    warning("Unknown metric. Using kendall.")
    metric <- "kendall"
  }

  worker <- function(z)
    dq.vect(z, min.cor, min.dist, q, metric, allgenes, min.genes,
            stratify.ds = stratify.ds, weighted.genes = weighted.genes,
            verbose = verbose)

  ## Single-vector shortcut
  if (is.numeric(ivect) && is.null(dim(ivect))) return(worker(ivect))

  if (jobs <= 1) {
    retl <- apply(ivect, 2, worker)
  } else {
    if (jobs > parallel::detectCores()) jobs <- parallel::detectCores()
    retl <- parallel::mclapply(seq_len(ncol(ivect)),
                               function(j) worker(ivect[, j]),
                               mc.cores = jobs)
  }
  do.call(rbind, retl)
}


#' Classify bulk gene-expression samples using the K-Nearest Neighbours (KNN) method
#'
#' Estimates the most probable iCMS subtype for each sample by majority vote among the
#' nn most similar prototype columns, with votes weighted by similarity score.
#'
#' @param ivect Gene-expression matrix with gene symbols as row names and samples as columns.
#'   Values should be log-transformed.
#' @param nn Number of nearest neighbours (default 10).
#' @param metric Similarity metric: "pearson", "kendall" (default), "spearman", or "cosine".
#' @param allgenes Use extended gene set (default FALSE; uses epithelial genes only).
#' @param min.genes Minimum gene overlap required (default 30).
#' @param mcenter Mean-centre each sample before computing cosine similarity (default FALSE;
#'   has no effect for correlation-based metrics).
#' @param stratify.ds If TRUE, vote is stratified by training dataset, preventing a single
#'   batch-matched cohort from dominating the top-nn neighbours (default FALSE).
#' @param weighted.genes If TRUE, up-weight discriminative genes by their absolute
#'   t-statistic when metric is "pearson" or "cosine"; for rank-based metrics restricts
#'   to the top 100 most discriminative genes (default FALSE).
#' @param jobs Number of parallel cores. Set to 1 to disable parallelism (default 4).
#' @param verbose If TRUE, include individual neighbour labels, similarities, and dataset
#'   labels in the output (default FALSE).
#' @return A data frame with one row per sample. Key columns: nearest.icms, confident.icms.
#' @export
iCMS.KNN <- function(ivect,
                     nn            = 10,
                     metric        = "kendall",
                     allgenes      = FALSE,
                     min.genes     = 30,
                     mcenter       = FALSE,
                     stratify.ds   = FALSE,
                     weighted.genes = FALSE,
                     jobs          = 4,
                     verbose       = FALSE) {

  message("input matrix: ", ncol(ivect), " samples, ", nrow(ivect),
          " genes. KNN method, metric=", metric,
          if (stratify.ds) ", stratify.ds=TRUE" else "",
          if (weighted.genes) ", weighted.genes=TRUE" else "")

  valid.metrics <- c("pearson", "kendall", "spearman", "cosine")
  if (!metric %in% valid.metrics) {
    warning("Unknown metric. Using kendall.")
    metric <- "kendall"
  }

  worker <- function(z)
    knn.vect(z, nn, metric, allgenes, min.genes,
             mcenter = mcenter, stratify.ds = stratify.ds,
             weighted.genes = weighted.genes, verbose = verbose)

  ## Single-vector shortcut
  if (is.numeric(ivect) && is.null(dim(ivect))) return(worker(ivect))

  if (jobs <= 1) {
    retl <- apply(ivect, 2, worker)
  } else {
    if (jobs > parallel::detectCores()) jobs <- parallel::detectCores()
    retl <- parallel::mclapply(seq_len(ncol(ivect)),
                               function(j) worker(ivect[, j]),
                               mc.cores = jobs)
  }
  do.call(rbind, retl)
}


#' Classify bulk gene-expression samples using the Nearest Template Prediction (NTP) method
#'
#' This is the published bulk iCMS classifier from Joanito et al. (Nat Genet 2022).
#' Note: NTP is sensitive to batch effects; consider batch correction before use.
#'
#' @param ivect Gene-expression matrix with gene symbols as row names and samples as columns.
#'   Values should be log-transformed.
#' @param metric Similarity metric (default "kendall").
#' @param bs.iter Number of bootstrap iterations for p-value estimation (default 100).
#' @param jobs Number of parallel cores (default 8).
#' @return A data frame with per-sample iCMS calls and bootstrap p-values.
#' @export
iCMS.NTP <- function(ivect, metric = "kendall", bs.iter = 100, jobs = 8) {
  ntp.matrix(ivect, tmat, metric, bs.iter, jobs)
}
