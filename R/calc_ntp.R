ntp.matrix <- function(imat, tmat, metric="kendall", bs.iter=100, jobs=8){

  if (jobs>detectCores()) { jobs <- detectCores() }
  cg <- intersect(rownames(tmat), rownames(imat))
  tmat <- tmat[cg,]

  isim <- apply(imat[cg,],2,function(x) vector.similarity(x,tmat, metric))
  isim <- t(isim)
  colnames(isim) <- c("i2","i3")

  ## bootstrap estimate of distance for random genes
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

  tmp <- unlist(mclapply(rownames(isim),
                         function(x) bootstrap.vector(imat[,x], tmat, isim[x,], metric, bs.iter)))
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
