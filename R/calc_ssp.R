#' calculate iCMS on bulk data (microarray / RNA-seq) using single-sample predictor
#' @param v input vector to calculate.
#' @param min.cor the minimum correlation threshold to determine the confidence
#' @param min.dist the minimum distance to determine the iCMS call
#' @param q The quantile that is used to judge distance. q=0.9-1 appears better for non batch and corrected data while q=0.6-0.7 appears better for batch-corrected
#' @param matric The method to calculate correlation matrix. Kendall correlation (non parametric) seems more robust but is slow. Pearson is faster, and often good enough.
#' @param allgenes  whether to use the extended genes or not, basically only epithelial genes are used
#' @param min.genes Minimum number of genes that can calculate iCMS call
#' @return iCMS classification for this sample
#' @export
#'
dq.vect <- function(v, min.cor, min.dist, q,  metric, allgenes, min.genes){
  if (allgenes) {
    prot <- ssproto.ext } else { prot <- ssproto }
  cg <- intersect(rownames(prot), names(na.omit(v)))
  if (length(cg)<min.genes) { warning("Too few common genes. Unable to map prototype") ; return(NA) }

  if (metric!="cosine"){
    smat <- cor(as.matrix(v)[cg,], prot[cg,], method=metric)
    colnames(smat) <- colnames(prot)
  } else {
    if (metric=="cosine") {
      cosine.sim <- function(x, y) {
        a <- crossprod(x,y)
        b <- outer(sqrt(apply(x,2,crossprod)), sqrt(apply(y,2,crossprod)))
        a/b
      }
      smat <- cosine.sim(as.matrix(v)[cg,], prot[cg,])
      colnames(smat) <- colnames(prot)
    }
  }

  maxdist <- t(apply(smat, 1, function(x) tapply(x, ssicms, quantile, q, na.rm=TRUE)))
  ## less robust?
  ## maxdist <- t(apply(smat, 1, function(x) tapply(x, ssicms, max, na.rm=TRUE)))
  probcms <- cbind(data.frame(smat),maxdist)
  probcms$i2i3 <- probcms$i2-probcms$i3

  ## nested function
  cfun <- function(x) {
    if (all(is.na(x))) { return(NA) }
    ## negative correlation with both prototypes
    if (max(x, na.rm=TRUE)<0) { return(NA) }
    ## cor must be > min.cor OR
    ## abs difference in cor between i2.i3
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

  return(probcms)
}

#' calculate cosine similarity
#' @export
cosine.sim <- function(x, y) {
  a <- crossprod(x,y)
  b <- outer(sqrt(apply(x,2,crossprod)), sqrt(apply(y,2,crossprod)))

  return(a/b)
}


