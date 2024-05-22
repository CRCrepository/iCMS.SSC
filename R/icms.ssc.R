#' calculate iCMS on bulk data (microarray / RNA-seq)
#' @param ivect input matrix to calculate. The row will be the gene symbol, and column will be sample name
#' @param min.cor the minimum correlation threshold to determine the confidence (default = 0.1)
#' @param min.dist the minimum distance to determine the iCMS call (default = 0.05)
#' @param q The quantile that is used to judge distance
#' @param matric The method to calculate correlation matrix (default = kendall)
#' @param allgenes description
#' @param min.genes Minimum number of genes that can calculate iCMS call (default = 30 genes)
#' @param jobs
#' @return iCMS classification
#' @export
iCMS.SSC <- function(ivect, min.cor=0.1, min.dist=0.05, q=0.9,
                     metric="kendall", allgenes=FALSE, min.genes=30, jobs=8) {

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

  if (class(ivect)[1]=="numeric") { return(ssp.vect(ivect, min.cor, min.dist,
                                                    q, metric, allgenes,
                                                    min.genes)) }

  if (jobs<=1) {
    retl <- apply(ivect, 2, ssp.vect, min.cor, min.dist, q,
                  metric, allgenes, min.genes)
  } else {
    if (jobs>detectCores()) { jobs <- detectCores() }
    retl <- mclapply(1:ncol(ivect), function(z) ssp.vect(ivect[,z], min.cor, min.dist, q,
                                                         metric, allgenes, min.genes))
  }

  do.call(rbind, retl)
}




ssp.vect <- function(v, min.cor, min.dist, q, metric=metric, allgenes, min.genes)
{
  if (allgenes) {
    prot <- ssproto.ext } else { prot <- ssproto }
  cg <- intersect(rownames(prot), names(na.omit(v)))
  if (length(cg)<min.genes) { warning("Too few common genes. Unable to map prototype") ; return(NA) }

  if (metric!="cosine"){
    smat <- cor(as.matrix(v)[cg,], prot[cg,], method=metric)
    colnames(smat) <- colnames(prot)
  } else {
    if (metric=="cosine") {
      smat <- cosine.sim(as.matrix(v)[cg,], prot[cg,])
      colnames(smat) <- colnames(prot)
    }
  }

  maxdist <- t(apply(smat, 1, function(x) tapply(x, ssicms.index, quantile, q, na.rm=TRUE)))
  ## less robust?
  ## maxdist <- t(apply(smat, 1, function(x) tapply(x, ssicms, max, na.rm=TRUE)))
  probcms <- cbind(data.frame(smat),maxdist)
  probcms$i2i3 <- probcms$i2-probcms$i3

  cfun <- function(x) {
    if (all(is.na(x))) { return(NA) }
    ## negative correlation with both prototypes
    if (max(x, na.rm=TRUE)<0) { return(NA) }
    ## cor must be > min.cor OR
    ## abs difference in cor betwen i2.i3
    ## must be > min.cor to call confidently
    if (max(x,na.rm=TRUE)<min.cor &
        (abs(x["i2"]-x["i3"])<min.cor))
    { return(NA) }
    ## dist from diag
    dd <- sqrt(2*(x["i2"]-x["i3"])^2)
    if (dd<min.dist) { return(NA) }
    colnames(maxdist)[which.max(x)]
  }

  probcms$nearest.icms <- apply(maxdist, 1, function(x) { if (all(is.na(x))) { return(NA) } ; colnames(maxdist)[which.max(x)] } )
  probcms$confident.icms <- apply(maxdist,1,cfun)
  probcms$nearest.dset <- apply(smat, 1, function(x) { if (all(is.na(x))) { return(NA) } ; colnames(smat)[which.max(x)] } )
  probcms$nearest.dset <- gsub("-s[0-9]+$","",probcms$nearest.dset)
  probcms$nearest.dset <- gsub("i[0-9]+[.]","",probcms$nearest.dset)
  probcms$ngenes <- length(cg)
  probcms
}


cosine.sim <- function(x, y) {
  a <- crossprod(x,y)
  b <- outer(sqrt(apply(x,2,crossprod)), sqrt(apply(y,2,crossprod)))
  a/b
}
