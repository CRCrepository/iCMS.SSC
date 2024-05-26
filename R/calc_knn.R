#' calculate iCMS on bulk data (microarray / RNA-seq) using k-nearest neighbors algorithm for each sample
#' @param v input matrix to calculate. The row will be the gene symbol, and column will be sample name
#' @param nn the minimum correlation threshold to determine the confidence (default = 0.1)
#' @param matric The method to calculate correlation matrix (default = kendall). Kendall correlation (non parametric) seems more robust but is slow. Pearson is faster, and often good enough.
#' @param allgenes  whether to use the extended genes or not, basically only epithelial genes are used (default = FALSE)
#' @param min.genes Minimum number of genes that can calculate iCMS call (default = 30 genes)
#' @param mcenter description
#' @return iCMS classification
#' @export
#'
knn.vect <- function(v, nn, metric, allgenes, min.genes, mcenter=FALSE, verbose=FALSE)
{
  if (allgenes) {
    prot <- ssproto.ext } else { prot <- ssproto }
  cg <- intersect(rownames(prot), names(na.omit(v)))

  if (length(cg)<min.genes) {
    warning("Too few common genes. Unable to map prototype")
    return(NA)
  }

  if (mcenter) {
    v <- v-mean(v, na.rm=TRUE)
  }

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

  ord <- order(smat, decreasing=TRUE)[1:nn]
  knn <- ssicms[ord]
  knnd <- as.vector(smat[ord])
  names(knnd) <- colnames(smat)[ord]
  knndset <- gsub("^i[23][.]|-s[0-9]+$","",colnames(smat)[ord])

  qd <- tapply(knnd, knn, median)
  i2 <- sum(knn=="i2")/nn
  i3 <- sum(knn=="i3")/nn
  i2i3 <- i2-i3

  cfun <- function(i2, i3) {
    if (i2>0.9) { return("i2") }
    if (i3>0.9) { return("i3") }
    return(NA)
  }

  nfun <- function(i2, i3, qd) {
    if (i2==i3) {
      if (qd["i2"]>qd["i3"]) { return("i2") }
      return("i3")
    }
    if (i2>i3) { return("i2") }
    return("i3")
  }


  if (verbose) {
    retl <- list(nn=knn,
                 nnd=knnd,
                 nnds=knndset)
  } else { retl <- list() }

  retl <- c(retl, qi2=as.vector(qd["i2"]),
            qi3=as.vector(qd["i3"]),
            i2=i2,
            i3=i3,
            i2i3=i2i3,
            nearest.icms=nfun(i2, i3, qd),
            confident.icms=cfun(i2,i3),
            nearest.ds=knndset[1],
            ngenes=length(cg))

  return(data.frame(retl))
}
