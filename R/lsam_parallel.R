library('clue')
library('base')
lsam_parallel <- function(mat){
  
  
}

lsam_classical <- function(mat){
  n <- nrow(mat)
  lsam.env <- new.env()
  assign('assign_row', rep(-1, n), lsam.env)
  assign('assign_col', rep(-1, n), lsam.env)
  # row reduction
  assign('dual_row', c(), lsam.env)
  for(i in seq(1, n)){
    lsam.env$dual_row <- c(lsam.env$dual_row, min(mat[i, ]))
  }
  assign('dual_col', c(), lsam.env)
  # column reduction
  for(j in seq(1, n)){
    lsam.env$dual_col <- c(lsam.env$dual_col, min(mat[,j] - lsam.env$dual_row))
  }
  check <- function(){
    match_count = 0 # check
    assign('cover_row', rep(0, n), lsam.env)
    assign('cover_col', rep(0, n), lsam.env)
    assign('prev_row', rep(-1, n), lsam.env)
    assign('prev_col', rep(-1, n), lsam.env)
    for(i in seq(1, n)){
      if(lsam.env$assign_row[i] != -1){
        lsam.env$cover_row[i] = 1
        match_count = match_count + 1
      }
    }
    if(match_count == n){
      #print(lsam.env$assign_col)
      return(lsam.env$assign_row) #Exit
    }
    st <- c() # Empty stack
    while(TRUE){ # [** search ** ]
      adj_list <- vector(mode='list',length=n) # initialize adjacency list
      for(i in seq(1, n)){
        if(lsam.env$cover_row[i] == 0)
          st <- c(st, i)
        for(j in seq(1, n)){
          if(mat[i, j] - lsam.env$dual_row[i] - lsam.env$dual_col[j])
            adj_list[[i]] <- c(adj_list[[i]], j)
        }
      }
      while(length(st)){ 
        i <- tail(st, n=1) # st.top
        st <- head(st, -1) # st.pop
        while(length(adj_list[[i]])){
          j <- adj_list[[i]][1]
          adj_list[[i]] <- tail(adj_list[[i]], -1)
          i_new <- lsam.env$assign_col[j]
          if(i_new == i)
            next()
          if(lsam.env$cover_col[j] == 0){
            lsam.env$prev_col[j] = i
            if(i_new == -1){
              # Augment
              # Procedure for augmenting current assignments by 1
              c_cur <- j
              r_cur <- -1
              while(c_cur != -1){
                r_cur <- lsam.env$prev_col[c_cur]
                lsam.env$assign_row[r_cur] <- c_cur
                lsam.env$assign_col[c_cur] <- r_cur
                c_cur <- lsam.env$prev_row[r_cur]
              }
              # Augment End
              return(check())
            }
            else{
              st <- c(st, i_new)
              lsam.env$prev_row[i_new] <- j
              lsam.env$cover_row[i_new] <- 0
              lsam.env$cover_col[j] <- 1
            }
          }
        }
      }
      # Update
      # Procedure for updating the dual variables
      theta <- Inf
      for(i in seq(1, n)){
        if(lsam.env$cover_row[i] == 0){
          for(j in seq(1, n)){
            if(lsam.env$cover_col[j] == 0){
              theta <- min(theta, mat[i, j] - lsam.env$dual_row[i] - lsam.env$dual_col[j])
            }
          }
        }
      }
      for(k in seq(1, n)){
        lsam.env$dual_row[k] <- lsam.env$dual_row[k] + theta / 2 * ifelse(lsam.env$cover_row[k] == 0, 1, -1)
        lsam.env$dual_col[k] <- lsam.env$dual_col[k] + theta / 2 * ifelse(lsam.env$cover_col[k] == 0, 1, -1)
      }
      # Update end
      # goto search
    }
  }
  return(check())
}

lsam_tree <- function(mat){
  n <- nrow(mat)
  lsam.env <- new.env()
  assign('assign_row', rep(-1, n), lsam.env)
  assign('assign_col', rep(-1, n), lsam.env)
  # row reduction
  assign('dual_row', c(), lsam.env)
  for(i in seq(1, n)){
    lsam.env$dual_row <- c(lsam.env$dual_row, min(mat[i, ]))
  }
  assign('dual_col', rep(Inf, n), lsam.env)
  # column reduction
  for(j in seq(1, n)){
    for(i in seq(1, n)){
      lsam.env$dual_col[j] <- min(lsam.env$dual_col, mat[i, j] - lsam.env$dual_row[i])
    }
  }
  check <- function(){
    match_count = 0 # check
    assign('cover_row', rep(0, n), lsam.env)
    assign('cover_col', rep(0, n), lsam.env)
    assign('prev_row', rep(-1, n), lsam.env)
    assign('prev_col', rep(-1, n), lsam.env)
    for(i in seq(1, n)){
      if(lsam.env$assign_row[i] != -1){
        lsam.env$cover_row[i] = 1
        match_count = match_count + 1
      }
    }
    if(match_count == n)
      return(lsam.env$assign_row) #Exit
    slack <- c()
    for(j in seq(1, n)){
      slack <- c(slack, Inf)
    }
    st <- c() # Empty stack
    for(i in seq(1, n)){
      if(lsam.env$cover_row[i] == 0)
        st <- c(st, i)
    }
    while(TRUE){ # [** search ** ]
      while(length(st)){ 
        i <- tail(st, n=1) # st.top
        st <- head(st, -1) # st.pop
        for(j in seq(1, n)){
          if(slack[j] > mat[i, j] - lsam.env$dual_row[i] - lsam.env$dual_col[j]){
            slack[j] <- mat[i, j] - lsam.env$dual_row[i] - lsam.env$dual_col[j]
            lsam.env$prev_col[j] <- i
          }
          if(!(mat[i, j] - lsam.env$dual_row[i] - lsam.env$dual_col[j])){
            i_new <- lsam.env$assign_col[j]
            if(lsam.env$cover_col[j] == 0){
              if(i_new == -1){
                # Augment
                # Procedure for augmenting current assignments by 1
                c_cur <- j
                r_cur <- -1
                while(c_cur != -1){
                  r_cur <- lsam.env$prev_col[c_cur]
                  lsam.env$assign_row[r_cur] <- c_cur
                  lsam.env$assign_col[c_cur] <- r_cur
                  c_cur <- lsam.env$prev_row[r_cur]
                }
                # Augment End
                return(check())
              }
              else{
                st <- c(st, i_new)
                lsam.env$prev_row[i_new] <- j
                lsam.env$cover_row[i_new] <- 0
                lsam.env$cover_col[j] <- 1
              }
            }
          }
        }
      }
      # Update2
      # Procedure for updating the dual soluation
      theta <- Inf
      for(j in seq(1, n)){
        theta <- min(theta, slack[j])
      }
      for(k in seq(1, n)){
        lsam.env$dual_row[k] <- lsam.env$dual_row[k] + theta / 2 * ifelse(lsam.env$cover_row[k] == 0, 1, -1)
        lsam.env$dual_col[k] <- lsam.env$dual_col[k] + theta / 2 * ifelse(lsam.env$cover_col[k] == 0, 1, -1)
      }
      for(j in seq(1, n)){
        if(slack[j] > 0){
          slack[j] <- slack[j] - theta
          if(slack[j] == 0){
            st <- c(st, lsam.env$prev_col[j])
          }
        }
      }
      # Update end
      # goto search
    }
  }
  return(check())
}




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
  return(mat)
}

"
mat=cost_mat(rnorm(100) * 100, rnorm(100) * 100)

start.time <- Sys.time()
a=lsam_classical(mat)
print(a)
end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)

start.time <- Sys.time()
b=lsam_tree(mat)
print(b)
end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)

start.time <- Sys.time()
c=solve_LSAP(mat)
print(c)
end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)
"

#mat=cost_mat(c(1,2,3), c(5,9,2))
#mat<-readRDS("mat-25.rds")
#print(mat)
#print(solve_LSAP(mat))
#print(lsam_classical(mat))
#print(lsam_tree(mat))
dim_mat<-1000
mat<-generate_mat(dim_mat)
start.time <- Sys.time()
c=solve_LSAP(mat)
print(c)
end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)

start.time <- Sys.time()
write.table(mat, file = "test.txt", sep = " ", row.names=FALSE, col.names=FALSE)
system('echo 1000 "$(cat test.txt)" | cpp_version/./hungarian', intern=TRUE)
end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)

