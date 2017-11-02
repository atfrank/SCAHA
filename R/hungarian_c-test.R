library('clue')
source('hungarian.R')

cost_mat <- function(x, y){
  xmat <- matrix(x,nrow=length(y),ncol=length(x),byrow=T)
  ymat <- matrix(y,nrow=length(y),ncol=length(x),byrow=F)
  costmat <- abs(xmat-ymat)
  return(costmat)
}
generate_mat <- function(n){
  mat=cost_mat(rnorm(n) * 100, rnorm(n) * 100)
  #mat[] <- vapply(mat, function(x){round(x)}, numeric(1))
  write.table(mat, file = paste("mat",n,".txt",sep=""), sep = " ", row.names = FALSE, col.names = FALSE)
  return(as.matrix(mat))
}

test_hungarian <- function(){
  dim_mat<-500
  mat<-generate_mat(dim_mat) # strange bug, works if write file and used read file, but not from direct matrix
  write.table(mat, file="test.txt", row.names = FALSE, col.names = FALSE)
  #mat<-as.matrix(read.table(file="test.txt"))
  colnames(mat) <- NULL
  get_cost <- function(xy){
    sum = 0
    for(i in 1:length(xy))
      sum = sum + mat[i, xy[i]]
    return(sum)
  }
  
  start.time <- Sys.time()
  a<-hungarian(mat)
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  
  start.time <- Sys.time()
  b<-solve_LSAP(mat)
  end.time <- Sys.time()
  time1.taken <- end.time - start.time
  
  #print(a)
  print("Custom Hungarian Algorithm:")
  print(time.taken)
  print(get_cost(a))
  #print(b)
  print("solve_LSAP:")
  print(time1.taken)
  print(get_cost(b))
}

hungarian <-
  function(x, maximum = FALSE)
  {
    if(!is.matrix(x) || any(x < 0))
      stop("x must be a matrix with nonnegative entries.")
    
    nr <- nrow(x)
    nc <- ncol(x)
    if(nr > nc)
      stop("x must not have more rows than columns.")
    if(nc > nr)
      x <- rbind(x, matrix(2 * sum(x), nc - nr, nc))
    
    if(maximum) x <- max(x) - x
    
    #TODO better support of double
    # idea: multiply all values by 1000 (3 decimal places), go back to int
    out <- .C("hungarian", n=as.integer(nc), x=as.integer(round(c(t(x)) * 1000)), p=integer(nc))$p + 1
    out <- out[seq_len(nr)]
    class(out) <- "hungarian"
    out
  }

print.hungarian <-
  function(x, ...)
  {
    writeLines(c("Optimal assignment:",
                 gsub("x", " ",
                      strwrap(paste(seq_along(x), x,
                                    sep = "x=>x", collapse = ", ")))))
    invisible(x)
  }

test_hungarian()
