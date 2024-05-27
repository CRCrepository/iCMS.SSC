#' applying NTP algorithm
#' @export
#'
ntp.matrix <- function(ivect, tmat, metric="kendall", bs.iter=100, jobs=4){

  if (jobs>parallel::detectCores()) { jobs <- parallel::detectCores() }
  cg <- intersect(rownames(tmat), rownames(ivect))
  tmat <- tmat[cg,]

  isim <- apply(ivect[cg,],2,function(x) vector.similarity(x,tmat, metric))
  isim <- t(isim)
  colnames(isim) <- c("i2","i3")

  ## bootstrap estivecte of distance for random genes
  bootstrap.vector <- function(v, tmat, simrow, metric=metric, bs.iter) {
    gc()
    rmat <- do.call(cbind,lapply(1:bs.iter, function(x) {
      v2 <- na.omit(v)
      if (length(v2)<nrow(tmat)) { return(NA) }
      v2[sample(1:length(v2), nrow(tmat))] } ))
    colnames(rmat) <- paste0("rsample",1:bs.iter)
    bootstrap.sim <- matrix.similarity(rmat, tmat, metric)
    sum(bootstrap.sim[,which.max(simrow)]>max(simrow))/bs.iter
  }

  if(jobs > 1){tmp <- unlist(parallel::mclapply(rownames(isim), function(x) bootstrap.vector(ivect[,x], tmat, isim[x,], metric, bs.iter), mc.cores = jobs))
  }else{
    tmp <- as.numeric( sapply(rownames(isim),function(x) bootstrap.vector(ivect[,x], tmat, isim[x,], metric, bs.iter)))
  }
  isim <- data.frame(isim)
  isim$nearest.icms <- apply(isim, 1, function(x) colnames(isim)[which.max(x)])
  isim$p.value <- tmp
  isim$confident.icms <- ifelse(isim$p.value<0.05, isim$nearest.icms, NA)
  isim$i2i3 <- isim$i2-isim$i3
  isim$icms.mixed <- isim$confident.icms
  isim$icms.mixed[is.na(isim$icms.mixed)] <- "i2.i3"
  isim$icms.max <- apply(isim[,c("i2","i3")], 1, max, na.rm=TRUE)

  return(isim)
}


#' calculate vector similarity #
#' @export
#'
vector.similarity <- function(v, tmat, metric)
{
  if (metric=="correlation") { metric  <- "pearson" }
  if (length(v)!=nrow(tmat)) {
    warning("Must have same number of genes")
    return(NaN)
  }

  cosine.sim <- function(x, y) {
    a <- crossprod(x,y)
    b <- outer(sqrt(apply(as.matrix(x),2,crossprod)), sqrt(apply(y,2,crossprod)))
    a/b
  }

  tmat <- tmat[!is.na(v),]
  v <- v[!is.na(v)]

  if (metric %in% c("pearson","kendall","spearman")) {
    smat <- cor(v, tmat, method=metric)
    colnames(smat) <- colnames(tmat)
  } else {
    if (metric=="cosine") {
      smat <- cosine.sim(v, tmat)
      colnames(smat) <- colnames(tmat)
    }
  }

  return(smat)
}


#' calculate matrix similarity
#' @export
matrix.similarity <- function(rmat, tmat, metric="kendall")
{
  if (metric=="correlation") { metric  <- "pearson" }
  if (nrow(rmat)!=nrow(tmat)) {
    warning("Must have same number of genes")
    return(NaN)
  }

  cosine.sim <- function(x, y) {
    a <- crossprod(x,y)
    b <- outer(sqrt(apply((x),2,crossprod)), sqrt(apply(y,2,crossprod)))
    a/b
  }

  if (metric %in% c("pearson","kendall","spearman")){
    smat <- cor(rmat, tmat, method=metric)
    colnames(smat) <- colnames(tmat)
  } else {
    if (metric=="cosine") {
      smat <- cosine.sim(rmat, tmat)
      colnames(smat) <- colnames(tmat)
    }
  }

  return(smat)
}

