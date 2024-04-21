loadPeaksData <- function(name, peaks.file) {
  if(!exists(name)) {
    rda.file = paste(name,".rda",sep="")
    if (file_test("-f", rda.file)) {
      load(rda.file,envir=parent.frame())
    }
    else {
      f <- file(peaks.file)
      test <- scan(file = f, what = "character", nmax = 1)
      close(f)
      if(test == "title") {
        peaks <- read.delim(peaks.file)
        peaks <- peaks[,-1]
        colnames(peaks) <- gsub("X","V",colnames(peaks))
        data <- as.matrix(peaks)
        assign(name, data, envir=parent.frame())
        save(list=c(name), file=rda.file, compress=TRUE, envir=parent.frame(n=1))
      }
      else
        loadCSVData(name, peaks.file)
    }
  }
  return(invisible(get(name,envir=parent.frame())))
}

loadCSVData <- function(name, csv.file) {
  csv <- read.csv2(csv.file, sep=";",header=FALSE)
  data <- as.matrix(csv)
  assign(name, data, envir=parent.frame())
  save(list=c(name), file=rda.file, compress=TRUE, envir=parent.frame(n=1))
  return(invisible(get(name,envir=parent.frame())))
}

ngsReadBED <- function(bed.file) {
  bed <- read.table(file=bed.file, stringsAsFactors=FALSE)
  if(bed[1,1] == "chr") {
    rnames <- bed[1,]
    bed <- bed[2:nrow(bed),]
    colnames(bed) <- rnames
  }
  bed
}

ngsSortData <- function(...,src.bed, dst.bed) {

  ## check
  if(nrow(src.bed) != nrow(dst.bed))
    stop("BED files have not the same size..")

  ## init.
  arg.sets <- list(...)
  names.sets <- as.character(as.list(substitute(list(...)))[-1L])
  
  ## compute ordering vector
  src.order <- order(src.bed[,1],src.bed[,2],src.bed[,3],src.bed[,6])
  dst.order <- order(dst.bed[,1],dst.bed[,2],dst.bed[,3],dst.bed[,6])
  orders <- data.frame(src=src.order,dst=dst.order)
  new.order <- orders[order(orders[,2]),1]

  ## sort sets
  sets <- list()
  for(i in 1:length(arg.sets)) {
    s <- arg.sets[[i]]
    s <- s[new.order,]
    if(length(sets) == 0)
      sets <- list(s)
    else
      sets <- c(sets,list(s))
  }
  
  ## eop
  names(sets) <- names.sets
  return(sets)
}

ngsSortBED <- function(src.bed, dst.bed) {

  ## check
  if(nrow(src.bed) != nrow(dst.bed))
    stop("BED files have not the same size..")

  ## compute ordering vector
  src.order <- order(src.bed[,1],src.bed[,2],src.bed[,3],src.bed[,6])
  dst.order <- order(dst.bed[,1],dst.bed[,2],dst.bed[,3],dst.bed[,6])
  orders <- data.frame(src=src.order,dst=dst.order)
  new.order <- orders[order(orders[,2]),1]

  ## eop
  return(new.order)
}

ngsWriteBED <- function(bed, out.file, order=NULL) {
  if(colnames(bed)[[1]]=="chr")
    cnames <- TRUE
  else
    cnames <- FALSE
  if(!is.null(order))
    bed <- bed[order,]
  write.table(bed,file=out.file,quote=FALSE,append=FALSE,sep="\t", 
              row.names=FALSE, col.names=cnames)
}
