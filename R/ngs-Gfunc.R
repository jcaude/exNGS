gfunc <- function(x,a,b,c) {
  y <- (a/b)*exp(-(x-c)^2/(2*b^2))
  y
}

gfunc.est <- function(xy,verbose=FALSE) {
  
  ## -- estimate A parameter
  a <- max(xy[WMIN:WMAX,"Y"])
  
  ## -- estimate C parameter
  a.idx <- which(xy[WMIN:WMAX,"Y"] == a)
  c <- a.idx[[which.min(abs(a.idx-WMID))]]+(WMAX-WMIN)/2-1

  if(verbose)
    cat("> A.est =", a, "C.est=",c,"\n")
  
  ## -- run non-linear regression analysis
  gcoefs <- try( nls(Y ~ (a/b)*exp(-(X-c)^2/(2*b^2)),
                     data=data.frame(xy),
                     start=list(a=a, b=2,c=c),
                     weights=W,
                     control=nls.control(maxiter=100),
                     trace=verbose),
                silent=TRUE)
  if(inherits(gcoefs,"try-error")) {
    a.hat = 0
    b.hat = 0.1
    c.hat = max(xy[,"X"])
    error = 1
  }
  else {
    model <- summary(gcoefs)
    a.hat <- coef(model)["a","Estimate"]
    b.hat <- coef(model)["b","Estimate"]
    c.hat <- coef(model)["c","Estimate"]
    error = 0
  }

  ## -- eop
  return(c(a.hat=a.hat, b.hat=b.hat, c.hat=c.hat, error=error))
}

gfunc.plot <- function(xy,gcoefs) {
  
  ## -- init
  model <- summary(gcoefs)
  a.hat <- coef(model)["a","Estimate"]
  b.hat <- coef(model)["b","Estimate"]
  c.hat <- coef(model)["c","Estimate"]
  
  ## -- console
  print(model)
  
  ## -- plot
  plot(xy)
  lines(WRANGE, gfunc(WRANGE, a.hat, b.hat, c.hat), col="red")
}

ngsGAUSS <- function(...,range,c.hat.range=15) {

  ## init
  add.sets <- list(...)
  names.sets <- as.character(as.list(substitute(list(...)))[-1L])
  ref <- add.sets[[1]]
  if(missing(range))
    range <- 1:ncol(ref)

  ## compute Gausian function estimators on the reference set
  gcoefs <- sfApply(ref, 1, function(x) { xy = data.frame(WRANGE,x) ;
                                          colnames(xy) = c("X","Y");
                                          g = gfunc.est(xy);
                                          g
                                        }
                    )
  
  ## sort reference set
  gf.order <-  order((1-gcoefs["error",])*gcoefs["a.hat",]*gcoefs["b.hat",], decreasing=TRUE)
  ref <- ref[gf.order,]
  rownames(ref) <- gf.order
  gcoefs <- gcoefs[,gf.order]
  
  # filter reference set               
  ref <- ref[(gcoefs["error",] == 0) & (abs(gcoefs["c.hat",]-100) < c.hat.range),range]
  ref.order <- as.numeric(rownames(ref))
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
