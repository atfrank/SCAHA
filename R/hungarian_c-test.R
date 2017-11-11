library('clue')
source('hungarian-test.R')

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
  dim_mat<-1000
  mat<-generate_mat(dim_mat) 
  #write.table(mat, file="test.txt", row.names = FALSE, col.names = FALSE)
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

test_hungarian()
