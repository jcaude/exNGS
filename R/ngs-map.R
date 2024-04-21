library(gplots)
library(ggplot2)
library(reshape2)
library(splines)

.splitMat <- function(mat, nset=1) {
  
  sets.start <- seq(1,ncol(mat), ncol(mat)/nset)
  sets.stop <- sets.start+(ncol(mat)/nset)-1
  sets.range <- cbind(start=sets.start, stop=sets.stop)
  
  sets <- list()
  for(i in 1:nrow(sets.range)) {
    start <- sets.range[i,"start"]
    stop <- sets.range[i,"stop"]
    range <- seq(start,stop,1)
    s <- mat[,range]
    if(length(sets) == 0)
      sets <- list(s)
    else
      sets <- c(sets,list(s))
  }
  return(sets)
}

.normFactor <- function(..., method=c("sum","max", "hi", "count"),cut=0) {
  
  ## init
  mat <- list(...)
  m <- mat[[1]]
  if(is.list(m) & length(mat) == 1)
    mat <- m
  
  ## compute normalization criteria
  crit <- list()
  for(i in 1:length(mat)) {
    x <- mat[[i]]
    x <- unlist(as.list(x))
    x <- x[x>cut]
    if(method == "sum")
      crit[i] <- sum(x)
    else if(method == "max")
      crit[i] <- max(x)
    else if(method == "hi")
      crit[i] <- fivenum(x)[4]
    else if(method == "count") 
      crit[i] <- length(mat[[i]]) - length(x)
    else
      crit[i] <- 1
  }

  ## eop
  return(crit)
}

ngsNormFactor <- function(..., method=c("none","sum.max","sum.avg", "max", "avg", "avg.g1")) {
  
  ## init
  mat <- list(...)
  m <- mat[[1]]
  if(is.list(m) & length(mat) == 1)
    mat <- m

  ## compute normalization factors
  if(method=="sum.max") {
    crit <- .normFactor(...,method="sum")
    target.crit <- max(unlist(crit))
  }
  if(method=="sum.avg") {
    crit <- .normFactor(...,method="sum")
    target.crit <- sum(unlist(crit))/length(mat)
  }
  else if(method == "max") {
    crit <- .normFactor(...,method="max")
    target.crit <- max(unlist(crit))
  }
  else if(method == "avg") {
    crit <- unlist(.normFactor(...,method="sum")) /
      unlist(.normFactor(...,method="count"))
    target.crit <- max(crit)
  }
  else if(method == "avg.g1") {
    crit <- unlist(.normFactor(...,method="sum",cut=1)) /
      unlist(.normFactor(...,method="count",cut=1))
    target.crit <- max(crit)
  }
  else
    return(rep.int(1,length(mat)))
  
  ## eop
  norm.factors <- sapply(crit, function(x) {target.crit/x})
  return(norm.factors)
}

ngsMap <- function(file, ..., max.rows=4000, norm.factors=NULL, scale='none') {

  ## init.
  sets <- list(...)
  s <- sets[[1]]
  if(is.list(s) & (length(sets) == 1))
    sets <- s
  if(length(sets) > 1) {
    csep <- sapply(sets,ncol)
    for(i in 2:length(csep))
      csep[i] <- csep[i-1] + csep[i]
  }
  else
    csep=NULL

  ## init. output
  if(grepl("pdf$", file))
    pdf(file=file, width=8.3, height=11.7)
  else
    png(file=file, width=8.3, height=11.7, units="in", res=150)

  ## DISPLAY NORMALIZATION

  ## - sum.avg by defaults 
  if(is.null(norm.factors) | (length(norm.factors) != length(sets))) 
    norm.factors <- ngsNormFactor(sets, method="sum.avg")
  for(i in 1:length(sets))
    sets[[i]] <- sets[[i]]*norm.factors[i]

  ## - build density matrix
  mat.density <- as.matrix(data.frame(sets))

  ## - cutoff (not for row scale)
  if(scale == 'none') {
    hinge <- min(unlist(.normFactor(sets, method="hi")))
    mat.density[mat.density > 5*hinge] <- 5*hinge
  }
  
  ## SAMPLE MATRIX (if necessary)
  if(nrow(mat.density) > max.rows) {
    samp.factor = nrow(mat.density) %/% max.rows
    mat.density = mat.density[seq(1,nrow(mat.density), samp.factor),]
  }
  
  ## build heatmap
  if(scale == 'none')
    cpanel <- colorpanel(256,"#FFFFFF","#66AAAA","#005555")
  else
    cpanel <- colorpanel(256,"#FFFFFF","#FFCCFF","#880088")
  ## ??? useRaster=TRUE -> disable Anti-Aliasing (perf++)
  heatmap.2(mat.density, scale=scale,
            dendrogram="none", Rowv=FALSE, Colv=FALSE,
            col=cpanel,
            key=FALSE, keysize=1.0, symkey=FALSE, density.info='none',
            trace='none', 
            colsep=csep, sepwidth=0.2,
            labRow="", labCol="",
            lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lhei=c(0.2, 4, 0.1 ), lwid=c(0.3,4)
            )

    
  ## eop
  dev.off()
}


ngsDensity <- function(file, ...) {

   ## init.
  sets <- list(...)
  s <- sets[[1]]
  if(is.list(s) & (length(sets) == 1))
    sets <- s

  ## init. output
  if(grepl("pdf$", file))
    pdf(file=file, width=8.3, height=11.7)
  else
    png(file=file, width=8.3, height=11.7, units="in", res=150)

  ## compute distributions
  raw.sum <- sapply(sets, function(x) { apply(x, 1, sum) })
  data <- data.frame(X=1:nrow(raw.sum), raw.sum)
  data <- melt(data,id="X")
  
  ## plot 
  p= ggplot(data, aes(X,value)) + geom_point(alpha=0.1)
  p= p + coord_flip() + facet_wrap(~ variable) + ylim(0,4000)
  p= p + stat_smooth(method="lm",formula=y~ns(x,3), colour="red")
  print(p)
  
  ## eop
  dev.off()
}


ngsProfile <- function(..., file, method=c("sum","mean")) {

  ## init.
  sets <- list(...)
  s <- sets[[1]]
  if(is.list(s) & (length(sets) == 1))
    sets <- s

  ## init. output
  if(!missing(file)) {
    if(grepl("pdf$", file))
      pdf(file=file, width=11.7, height=8.3)
    else
      png(file=file, width=11.7, height=8.3, units="in", res=150)
  }

  ## compute profile
  profile.sum <- sapply(sets, function(x) { apply(x, 2, method) })
  data <- data.frame(X=1:nrow(profile.sum), profile.sum)
  data <- melt(data,id="X")
  colnames(data) <- gsub("variable", "datasets", colnames(data))

  ## plot
  p = ggplot(data, aes(X,value,group=datasets)) + geom_line(aes(colour=datasets))
  p = p + ylab(paste(method,"of reads"))
  p = p + xlab("Position") + xlim(0,nrow(profile.sum))
  print(p)
  
  ## eop
  if(!missing(file))
    dev.off()
}
