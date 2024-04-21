.ngs <- function(...,
                 criteria=c("sum","mean","median","IQR","wsum","sum_range",
                            "wsum2","xmode"),
                 range,crit.only=FALSE) {

  ## local function
  gfunc <- function(x,a,b,c) {
    y <- (a/b)*exp(-(x-c)^2/(2*b^2))
    y <- y+1
    y
  }

  xmode <- function(x,weights) {
    mode <- length(x) - which.max(x)
    mode
  }
  
  wsum <- function(x,weights) {
    sum <- sum(weights*x)
    sum
  }
  
  wsum2 <- function(x,weights) {
    sum <- sum(weights*x) - max(x)
    
    idx <- which(x==max(x))
    mpos <- idx[which.min(abs(idx-WMID))]
    
    sum <- sum + sum(max(weights)-weights[max(0,mpos-1):min(ncol(x),mpos+1)])*sign(mpos-20)
    sum
  }

  ## init
  add.sets <- list(...)
  names.sets <- as.character(as.list(substitute(list(...)))[-1L])
  ref <- add.sets[[1]]
  if(missing(range)) range <- 1:ncol(ref)
  WMID <- ncol(ref)/2
  W <- gfunc(1:ncol(ref),60,5,WMID)

  ## init (special case)
  if(criteria == "sum_range") {
    calc.range = range
    range = 1:ncol(ref)
  }
  
  ## sort reference set
  if(criteria == "wsum")
    crit.ref <- apply(ref,1,wsum,W)
  else if(criteria == "wsum2")
    crit.ref <- apply(ref,1,wsum2,W)
  else if(criteria == "sum_range")
    crit.ref <- apply(ref[,calc.range], 1, sum)
  else 
    crit.ref <- apply(ref,1,criteria)

  ## if we only want the criteria value..
  if(crit.only) return(crit.ref)
  
  crit.order <- order(crit.ref, decreasing=TRUE)
  ref <- ref[crit.order,range]
  rownames(ref) <- crit.order
  sets <- list(ref)

  ## sort additional sets
  if(length(add.sets) > 1) {
    for(i in 2:length(add.sets)) {
      s <- add.sets[[i]]
      s <- s[crit.order,range]
      rownames(s) <- crit.order
      sets <- c(sets,list(s))
    }
  }
      
  ## eop
  names(sets) <- names.sets
  return(sets)
}

ngsSUMRange <- function(...,range) {
  if(missing(range))
    stop("Sum Range method requires a range value")
  else
    return(.ngs(...,criteria="sum_range",range=range))
}

ngsSUM <- function(...,range) {
  if(missing(range))
    return(.ngs(...,criteria="sum"))
  else
    return(.ngs(...,criteria="sum",range=range))
}

ngsMEAN <- function(...,range) {
  if(missing(range))
    return(.ngs(...,criteria="mean"))
  else
    return(.ngs(...,criteria="mean",range=range))
}

ngsMEDIAN <- function(...,range) {
  if(missing(range))
    return(.ngs(...,criteria="median"))
  else
    return(.ngs(...,criteria="median",range=range))
}

ngsWSUM <- function(...,range) {
  if(missing(range))
    return(.ngs(...,criteria="wsum"))
  else
    return(.ngs(...,criteria="wsum",range=range))
}

ngsWSUM2 <- function(...,range) {
  if(missing(range))
    return(.ngs(...,criteria="wsum2"))
  else
    return(.ngs(...,criteria="wsum2",range=range))
}

ngsIQR <- function(...,range) {
  if(missing(range))
    return(.ngs(...,criteria="IQR"))
  else
    return(.ngs(...,criteria="IQR",range=range))
}

ngsXMODE <- function(...,range) {
  if(missing(range))
    return(.ngs(...,criteria="xmode"))
  else
    return(.ngs(...,criteria="xmode",range=range))
}

ngsRND <- function(...,range) {

  ## init
  add.sets <- list(...)
  names.sets <- as.character(as.list(substitute(list(...)))[-1L])
  ref <- add.sets[[1]]
  if(missing(range))
    range <- 1:ncol(ref)

  ## sample first set
  rnd.order <- sample(1:nrow(ref), nrow(ref))
  ref <- ref[rnd.order,range]
  rownames(ref) <- rnd.order
  sets <- list(ref)

  ## sort additional sets
  if(length(add.sets) > 1) {
    for(i in 2:length(add.sets)) {
      s <- add.sets[[i]]
      s <- s[rnd.order,range]
      rownames(s) <- rnd.order
      sets <- c(sets,list(s))
    }
  }
  
  ## eop
  names(sets) <- names.sets
  return(sets)
}

ngsNONE <- function(...,range) {

  ## init
  add.sets <- list(...)
  names.sets <- as.character(as.list(substitute(list(...)))[-1L])
  ref <- add.sets[[1]]
  if(missing(range))
    range <- 1:ncol(ref)

  ## sample first set
  ref.order <- 1:nrow(ref)
  ref <- ref[ref.order,range]
  rownames(ref) <- ref.order
  sets <- list(ref)

  ## sort additional sets
  if(length(add.sets) > 1) {
    for(i in 2:length(add.sets)) {
      s <- add.sets[[i]]
      s <- s[ref.order,range]
      rownames(s) <- ref.order
      sets <- c(sets,list(s))
    }
  }

  ## eop
  names(sets) <- names.sets
  return(sets)
}
