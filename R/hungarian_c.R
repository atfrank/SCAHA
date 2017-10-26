library('clue')
dyn.load("hungarian.so")

cost_mat <- function(x, y){
  xmat <- matrix(x,nrow=length(y),ncol=length(x),byrow=T)
  ymat <- matrix(y,nrow=length(y),ncol=length(x),byrow=F)
  costmat <- abs(xmat-ymat)
  return(costmat)
}
generate_mat <- function(n){
  mat=cost_mat(rnorm(n) * 100, rnorm(n) * 100)
  mat[] <- vapply(mat, function(x){round(x)}, numeric(1))
  write.table(mat, file = paste("mat",n,".txt",sep=""), sep = " ", row.names = FALSE, col.names = FALSE)
  return(as.matrix(mat))
}

dim_mat<-500
mat<-generate_mat(dim_mat) # strange bug, works if write file and used read file, but not from direct matrix
write.table(mat, file="test.txt", row.names = FALSE, col.names = FALSE)
mat<-as.matrix(read.table(file="test.txt"))
colnames(mat) <- NULL
get_cost <- function(xy){
  sum = 0
  for(i in 1:length(xy))
    sum = sum + mat[i, xy[i]]
  return(sum)
}

test <- c(t(mat))
print(mat[1,1:100])
print(test[1:100])

start.time <- Sys.time()
a<-.C("hungarian", n=as.integer(dim_mat), x=test, p=integer(dim_mat))$p + 1
end.time <- Sys.time()
time.taken <- end.time - start.time

start.time <- Sys.time()
b<-solve_LSAP(mat)
end.time <- Sys.time()
time1.taken <- end.time - start.time

#print(a)
print(time.taken)
print(get_cost(a))
#print(b)
print(time1.taken)
print(get_cost(b))
