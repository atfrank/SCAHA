lsam_parallel <- function(mat){
  
  
}

lsam_classical <- function(mat){
  n <- nrow(mat)
  assign_row <- rep(-1, n)
  assign_col <- rep(-1, n)
  # row reduction
  dual_row <- c()
  for(i in seq(1, n)){
    dual_row <- c(dual_row, min(mat[i, ]))
  }
  dual_col <- c()
  # column reduction
  for(j in seq(1, n)){
    dual_col <- c(dual_col, min(mat[,j] - dual_row))
  }
  print(dual_row)
  print(dual_col)
  while(TRUE){
    match_count = 0 # check
    cover_row <- rep(0, n)
    cover_col <- rep(0, n)
    prev_row <- rep(-1, n)
    prev_col <- rep(-1, n)
    for(i in seq(1, n)){
      if(assign_row[i] != -1){
        cover_row[i] = 1
        match_count = match_count + 1
      }
    }
    if(match_count == n)
      break
    st <- c()
    # stack
    for(i in seq(1, n)){
      if(cover_row[i] == 0)
        st <- c(stack, i)
      adj_list <- vector(mode='list',length=n)
      for(j in seq(1, n)){
        if(mat[i, j] - dual_row[j] == 0)
          adj_list[i] <- c(adj_list[i], j)
      }
    }
    while(!is.null(st)){ # search
      i <- tail(st, n=1) # st.top
      st <- head(st, -1) # st.pop
      while(!is.null(adj_list[i])){
        j <- adj_list[i][[0]]
        adj[[i]] <- NULL
        i_new <- assign_col[j]
        if(i_new == i)
          next()
        if(cover_col[j] == 0){
          prev_col[j] = i
          if(i_new == -1){
            # Augment
            # Procedure for augmenting current assignments by 1
            c_cur <- j
            r_cur <- -1
            while(c_cur != -1){
              r_cur <- prev_col[c_cur]
              assign_row[r_cur] <- c_cur
              assign_col[c_cur] <- r_cur
              c_cur <- prev_row[r_cur]
            }
            # Augment End
            # TODO goto check 
          }
          else{
            st <- c(st, i_new)
            prev_row[i_new] <- j
            cover_row <- 0
            cover_col <- 1
          }
        }
      }
    }
    # Update
    # Procedure for updating the dual variables
    theta <- Inf
    for(i in seq(1, n)){
      if(cover_row[i] == 0){
        for(j in seq(1, n)){
          if(cover_col == 0){
            theta <- min(theta, mat[i, j] - dual_row[i] - dual_col[j])
          }
        }
      }
      for(k in seq(1, n)){
        if(cover_row[k] == 0){
          dual_row[k] <- dual_row[k] + theta / 2 * ifelse(cover_row[k] == 0, 1, -1)
        }
        if(cover_col[k] == 0){
          dual_col[k] <- dual_col[k] + theta / 2 * ifelse(cover_col[k] == 0, 1, -1)
        }
      }
    }
    #Update end
    # TODO goto search
  }
  return(assign_row)
}

lsam_tree <- function(){
  n <- nrow(mat)
  # Initial reduction
  assign_row <- rep(-1, n)
  assign_col <- rep(-1, n)
  # row reduction
  dual_row <- c()
  for(i in seq(1, n)){
    dual_row <- c(dual_row, min(mat[i, ]))
  }
  dual_col <- c()
  # column reduction
  for(j in seq(1, n)){
    dual_col <- c(dual_col, min(mat[,j] - dual_row))
  }
  while(TRUE){
    # optimality check
    match_count = 0
    cover_row <- rep(0, n)
    cover_col <- rep(0, n)
    prev_row <- rep(-1, n)
    prev_col <- rep(-1, n)
    for(i in seq(1, n)){
      if(assign_row[i] != -1){
        cover_row[i] = 1
        match_count = match_count + 1
      }
    }
    if(match_count == n)
      break
    # Augmenting path search
    st <- c()
    for(i in seq(1, n)){
      if(cover_row == 0)
        st <- c(st, i)
    }
    while(!is.null(st)){
      # search
      i <- tail(st, n=1) # top
      st <- head(st, -1) # pop
      for(j in seq(1, n)){
        
      }
    }
  }
}




cost_mat <- function(x, y){
  xmat <- matrix(x,nrow=length(y),ncol=length(x),byrow=T)  
  ymat <- matrix(y,nrow=length(y),ncol=length(x),byrow=F)
  costmat <- abs(xmat-ymat)
  return(costmat)
}

mat=cost_mat(c(1,2,3), c(5,9,2))
#lsam_classical(mat)


