#' calculate iCMS on bulk data (microarray / RNA-seq) using single-sample predictor
#' @param ivect input matrix to calculate. The row will be the gene symbol, and column will be sample name
#' @param min.cor the minimum correlation threshold to determine the confidence (default = 0.1)
#' @param min.dist the minimum distance to determine the iCMS call (default = 0.05)
#' @param q The quantile that is used to judge distance (default = 0.9). q=0.9-1 appears better for non batch and corrected data while q=0.6-0.7 appears better for batch-corrected
#' @param matric The method to calculate correlation matrix (default = kendall). Kendall correlation (non parametric) seems more robust but is slow. Pearson is faster, and often good enough.
#' @param allgenes  whether to use the extended genes or not, basically only epithelial genes are used (default = FALSE)
#' @param min.genes Minimum number of genes that can calculate iCMS call (default = 30 genes)
#' @param jobs If parallel execution fails you can try jobs=1 which will use serial apply (default = 4)
#' @return iCMS classification
#' @export
#'
iCMS.DQ <- function(ivect, min.cor=0.1, min.dist=0.05, q=0.9,
                     metric="kendall", allgenes=FALSE, min.genes=30, jobs=4) {

  message("input matrix : ",  ncol(ivect), " samples with ", nrow(ivect) ," genes. Use Single-Sample Predictor with ",metric," method")

  ## q is the quantile that is used to judge distance
  ## q=1 means use the max, ie the most similar centroid
  if (q>1) { q <- 1 }
  if (q<0) { q <- 0 }
  if (q<0.5) { warning("Are you sure you want q<0.5? This will consider the distance of samples with similarity below the median.") }

  valid.metrics <- c("pearson", "kendall", "spearman", "cosine")
  if (!metric %in% valid.metrics) {
    warning(paste0("Only understand metrics: ", valid.metrics))
    warning("Using Kendall correlation, by default")
    metric <- "kendall"
  }

  if (class(ivect)[1]=="numeric") { return(dq.vect(ivect, min.cor, min.dist,
                                                    q, metric, allgenes,
                                                    min.genes)) }

  if (jobs<=1) {
    retl <- apply(ivect, 2, function(z) dq.vect(z, min.cor, min.dist, q, metric, allgenes, min.genes))
  } else {
    if (jobs>detectCores()) { jobs <- detectCores() }
    retl <- parallel::mclapply(1:ncol(ivect), function(z) dq.vect(ivect[,z], min.cor, min.dist, q,
                                                         metric, allgenes, min.genes), mc.cores = jobs)
  }

  return(do.call(rbind, retl))
}


#' calculate iCMS on bulk data (microarray / RNA-seq) using k-nearest neighbors algorithm
#' @param ivect input matrix to calculate. The row will be the gene symbol, and column will be sample name
#' @param nn the minimum correlation threshold to determine the confidence (default = 0.1)
#' @param matric The method to calculate correlation matrix (default = kendall). Kendall correlation (non parametric) seems more robust but is slow. Pearson is faster, and often good enough.
#' @param allgenes  whether to use the extended genes or not, basically only epithelial genes are used (default = FALSE)
#' @param min.genes Minimum number of genes that can calculate iCMS call (default = 30 genes)
#' @param jobs If parallel execution fails you can try jobs=1 which will use serial apply (default = 4)
#' @return iCMS classification
#' @export
#'
iCMS.KNN <- function(ivect, nn=10,  metric="kendall", allgenes=FALSE, min.genes=30, jobs=4, verbose = FALSE) {

  message("input matrix : ",  ncol(ivect), " samples with ", nrow(ivect),". Use K-Nearest Neighbors algorithm with ",metric," method")

  valid.metrics <- c("pearson", "kendall", "spearman", "cosine")
  if (!metric %in% valid.metrics) {
    warning(paste0("Only understand metrics: ", valid.metrics))
    warning("Using Kendall correlation, by default")
    metric <- "kendall"
  }

  if (class(ivect)[1]=="numeric") { return(knn.vect(ivect, nn,  metric, allgenes,  min.genes, verbose)) }

  if (jobs<=1) {
    retl <- apply(ivect, 2, knn.vect, min.cor, min.dist, q,metric, allgenes, min.genes, verbose)
  } else {
    if (jobs>detectCores()) { jobs <- detectCores() }
    retl <- parallel::mclapply(1:ncol(ivect), function(z) knn.vect(ivect[,z], nn, metric, allgenes, min.genes, verbose), mc.cores = jobs)
  }

  return(do.call(rbind, retl))
}

#' calculate iCMS on bulk data (microarray / RNA-seq) using k-nearest neighbors algorithm
#' @param ivect input matrix to calculate. The row will be the gene symbol, and column will be sample name
#' @param matric The method to calculate correlation matrix (default = kendall). Kendall correlation (non parametric) seems more robust but is slow. Pearson is faster, and often good enough.
#' @param bs.iter the number of iteration to calculate P-value
#' @param allgenes  whether to use the extended genes or not, basically only epithelial genes are used (default = FALSE)
#' @param min.genes Minimum number of genes that can calculate iCMS call (default = 30 genes)
#' @param jobs If parallel execution fails you can try jobs=1 which will use serial apply (default = 4)
#' @return iCMS classification
#' @export
#'
iCMS.NTP <- function(ivect, metric="kendall", bs.iter=100, jobs=8){
  return(ntp.matrix(ivect, tmat, metric, bs.iter, jobs))
}

