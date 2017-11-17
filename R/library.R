# Functions
check_and_install <- function(pkg){
  # check if package is already installed
  # if not install
  if(!(pkg %in% installed.packages())){install.packages(pkgs = pkg)}
}  

install_scaha_dependencies <- function(){
  # optparse: for writing neat Unix-like commandline R script
  # plyr: for ddply algorithm
  # clue: for solve_LSAP (Hungarian algorithm)
  # doParallel: for running in parallel
  pkgs <- c("optparse", "plyr", "clue", "doParallel")
  for (pkg in pkgs){
    check_and_install(pkg)
  }
}

load_assigned_computed_shifts <- function(assgined_computed_cs_filename){
  # loads assigned computed chemical shift data file
  assgined_computed_cs <- read.table(assgined_computed_cs_filename, header=FALSE)  
  # add appropiate colname names
  if (ncol(assgined_computed_cs==5)){names <- c("conformation", "resid", "resname", "nucleus", "cs")}
  if (ncol(assgined_computed_cs==6)){names <- c("conformation", "resid", "resname", "nucleus", "cs", "na")}       
  colnames(assgined_computed_cs) <- names

  # only keep necessary columns from
  assgined_computed_cs <- assgined_computed_cs[,names]
  assgined_computed_cs$error <- 1
  levels <- length(unique(assgined_computed_cs$conformation))
  assgined_computed_cs['assigned'] <- NA
  return(assgined_computed_cs)
}

load_unassigned_shifts <- function(unassgined_data_file, testing){
  # loads unassigned chemical shift data file
  unassgined_data <- read.table(unassgined_data_file, header=FALSE)
  # in testing mode expects a single column file
  # just in case user supplies truly unassigned peaks but sets testing to TRUE
  if(ncol(unassgined_data) == 1){testing <- FALSE}
  # add appropiate colname names
  if (ncol(unassgined_data) == 5){colnames(unassgined_data) <- c("resname","resid","nucleus","cs","dummy")}
  if (ncol(unassgined_data) ==  1){colnames(unassgined_data) <- c("cs")}
  return (unassgined_data)
}

duplicate_unassigned_data <- function(unassgined_data, assgined_computed_cs){
  # duplicate data if data in assigned computed chemical shift list is greater than unassgined_data
  nassgined_computed_cs <- nrow(subset(assgined_computed_cs, conformation==1))
  nrefcs <- nrow(unassgined_data)
  if (nassgined_computed_cs > nrefcs){
    for (i in 1:ceiling(nassgined_computed_cs/nrefcs)){
      unassgined_data <- rbind(unassgined_data, unassgined_data)
    }
  }
  return(unassgined_data)
}

assign <- function(x, y, z=NULL, custom = FALSE, twoD = FALSE){
  suppressPackageStartupMessages(require("clue"))
  # assigns elements in y to elements in x 
  # Input -- x (vector): e.g., chemical shifts 
  # Input -- y (vector): e.g., chemical shifts which are to be mapped to those in x (note that x can be greater than y)
  # Input -- z (vector): optional weights  
  if (!twoD){
    xmat <- matrix(x,nrow=length(y),ncol=length(x),byrow=T)  
    ymat <- matrix(y,nrow=length(y),ncol=length(x),byrow=F)
    costmat <- abs(xmat-ymat)
  } else {
    costmat <- matrix(0, nrow = nrow(unassgined_data), ncol = nrow(assgined_computed_cs))
    for (i in 1:nrow(unassgined_data)){
      d1 <- matrix(as.vector(unlist(unassgined_data[i, c("ppm1", "ppm2")])), ncol=2, nrow=nrow(assgined_computed_cs), byrow = TRUE)
      d2 <- matrix(as.vector(unlist(assgined_computed_cs[, c("peak_H", "peak_C")])), ncol=2, byrow = FALSE)
      w2 <- matrix(as.vector(unlist(assgined_computed_cs[, c("weight_H", "weight_C")])), ncol=2, byrow = FALSE)
      costmat[i, ] <- w2[,1]*abs(d1[,1]-d2[,1])+w2[,2]*abs(d1[,2]-d2[,2])
    }
  }
  
  # should I use a weighted cost matrix or not (not yet tested)
  if (is.null(z)){
		if(custom){
      # use custom Hungarian algorithm
      a <- hungarian(costmat)
    }
    else{
      # use the Hungarian algorithm implemented in solve_LSAP
      a <- solve_LSAP(costmat)
    }
  } else {
		zmat <- matrix(1/z,nrow=length(y),ncol=length(x),byrow=T)
		if(custom){
		  # use custom Hungarian algorithm
		  a <- hungarian(costmat * zmat)
		}
		else{
		  # use the Hungarian algorithm implemented in solve_LSAP
		  a <- solve_LSAP(costmat * zmat)
		}
  }
  
  # make list containing the raw assignments, assignment matrix, and cost matrix
  # and then return
  a_mat <- matrix(0,nrow=length(y),ncol=length(x),byrow=T)
  for (i in 1:length(y)){
    j <- as.vector(a)[i]
    a_mat[i,j] <-1
  }
  return(list(a=a,a_mat=a_mat,costmat=costmat))
}

get_assignment_cost <- function(a, costmat){
  # get cost of an assignment (useful for debugging)
  # Input -- a (vector): linear assignment object returned by solve_LSAP
  # Input -- costmat (matrix): cost matrix
  cost = as.integer(0)
  for(i in 1:length(a))
    cost = cost + costmat[i, a[i]]
  return(cost)
}

get_assignments <- function(assgined_computed_cs, unassgined_data=unassgined_data, testing=FALSE, output=output, weighted=FALSE, iter=1, freq_output=1, scale=0.0, custom = FALSE, twoD = FALSE){
  # Assign Peaks Using the Hungarian Algorithm 
  # Input  -- assgined_computed_cs (dataframe): assigned computed chemical shifts (must contain columns: c("resid", "resname", "nucleus", "cs", "error"))
  # Input  -- unassgined_data (dataframe): unassigned chemical shift peak (expect this data under column named c("cs"))
  # Option -- testing (logical): true is the unassgined_data frame does contain the correct assignment; useful testing; In that case unassgined_data should contain columns: c("resid", "resname", "nucleus", "cs")
  # Option -- output (character string): base name for output assignmen files
  # Option -- weighted (logical): do weighted assignment?
  # Option -- iter (integer): number of assignments to carry out (default is 1 (this is what we used in paper)
  # Option -- freq_output (integer): frequency with which to write assignments; only valid if iter > 1
  # Option -- scale (double): scale factor determining amount of noise to add to computed chemical shift; zero corresponds to no added noise (default is 0 (this is what we used in paper)

  # read in assignment computed chemical shift data
  assgined_computed_cs <- assgined_computed_cs[order(assgined_computed_cs$resid, assgined_computed_cs$nucleus),]  
  if(!twoD){
    unassgined_data_only <- unassgined_data[,"cs"]
  } else {
    unassgined_data_only <- unassgined_data[,c("ppm1", "ppm2")]
  }
  
  # initialize matrix for a probability matrix
  # not need in the current implementation of SCAHA since we've set iter
  # but may useful, in the future, if we do multiple assignments; in that case we can return the most probable assignment)
  prob <- matrix(0, nrow=nrow(assgined_computed_cs), ncol=nrow(unassgined_data))
  conformation <- unique(assgined_computed_cs$conformation)   
  
  # looping over multiple assignments
  # however, in the current version of SCAHA, iter = 1
  for(i in 1:iter){
    # do actual assignment
    if (weighted){      
    	a <- assign(unassgined_data_only, assgined_computed_cs$cs, z=assgined_computed_cs$error, custom = custom, twoD = twoD)
    } else {
    	a <- assign(unassgined_data_only, assgined_computed_cs$cs, custom = custom, twoD = twoD)
    }
    # store a probablity matrix
    prob <- prob + a$a_mat
    
    if (i%%freq_output==0){
			# get most probable assignment 
			assgined_computed_cs_prob <- prob / i
			a <- solve_LSAP(1 - assgined_computed_cs_prob)
			# map unassigned chemical shift peaks to the assigned computed chemical shifts based on the
			# assignments determined by the Hungarian algorithm
			assgined_computed_cs$assigned <- unassgined_data[a,"cs"]	
			
			if (testing){
				# for testing: merge pseudo-unassigned chemical shift data with the predicted assignments
				# allow us to easily assess assignment accuracy for cases where we do know the correct assignment
				if(!twoD){
				  tmp <- merge(unassgined_data, assgined_computed_cs, by=c("resid","nucleus"))
				  tmp <- tmp[order(tmp$resid, tmp$nucleus),]
				  tmp <- rename(tmp, c("resname.x" = "resname", "cs.x"="actual", "cs.y"="predicted"))
				  tmp <- unique(tmp[, c("conformation", "resid", "resname", "nucleus", "actual", "assigned", "predicted")])			
				} else {
				  tmp <- merge(unassgined_data, assgined_computed_cs, by=c("resid","nucleus"))
				  tmp <- tmp[order(tmp$resid, tmp$nucleus),]
				  tmp <- rename(tmp, c("resname.x" = "resname", "cs.x"="actual", "cs.y"="predicted"))
				  tmp <- unique(tmp[, c("conformation", "resid", "resname", "nucleus", "actual", "assigned", "predicted")])			
				}
			} else {
				if(!twoD){
				  tmp <- rename(assgined_computed_cs, c("cs"="predicted"))
				  tmp <- unique(tmp[, c("conformation", "resid", "resname", "nucleus", "assigned", "predicted")])			
				} else {
				  tmp <- rename(assgined_computed_cs, c("cs"="predicted"))
				  tmp <- unique(tmp[, c("conformation", "resid", "resname", "nucleus", "assigned", "predicted")])			
				}
			}
			# output assignments to a text file
			write.table(tmp, file = paste(output, "_", i, "_", conformation, ".txt", sep = ""), quote = FALSE, row.names = FALSE, col.names = FALSE)
			if (verbose){print(tmp)}			
		}    
  }
  return(TRUE)
}

SCAHA <- function(assgined_computed_cs_filename="data/predicted_CS_table_test_clean.txt", unassgined_data=unassgined_data, testing=FALSE, iter=1, freq_output = 1, scale=1, output=output, parallel = FALSE, weighted=FALSE, conformation_one = FALSE, custom = FALSE){
  suppressPackageStartupMessages(require("plyr"))
  # Input  -- assgined_computed_cs_filename (character string): path to assigned computed chemical shift data
  # Input  -- unassgined_data_file (character string): path to  unassigned chemical shift peak (expects this data under column named c("cs"))
  # Option -- testing (logical): true is the unassgined_data frame does contain the correct assignment; useful testing; In that case unassgined_data should contain columns: c("resid", "resname", "nucleus", "cs")
  # Option -- iter (integer): number of assignments to carry out (default is 1 (this is what we used in paper)
  # Option -- freq_output (integer): frequency with which to write assignments; only valid if iter > 1
  # Option -- scale (double): scale factor determining amount of noise to add to computed chemical shift; zero corresponds to no added noise (default is 0 (this is what we used in paper)
  # Option -- output (character string): base name for output assignmen files
  # Option -- parallel (logical): run in parallel mode?
  # Option -- weighted (logical): do weighted assignment?
  # Option -- conformation_one (logical): run of conformation one? Useful for debugging
  
  # loads assigned computed chemical shift data
  assgined_computed_cs <- load_assigned_computed_shifts(assgined_computed_cs_filename)
  
  # load unassigned chemical shift data
  unassgined_data <- load_unassigned_shifts(unassgined_data_file, testing)
  
  # get data for along one conformation if conformation_one == TRUE; useful for testing
  if (conformation_one){assgined_computed_cs <- subset(assgined_computed_cs, conformation==1)}
  
  # duplicate data if data in assigned computed chemical shift list is greater than unassgined_data
  unassgined_data <- duplicate_unassigned_data(unassgined_data, assgined_computed_cs)
  
  # run get_assignments (work-horse)
  ddply(.data=assgined_computed_cs, .var=c("conformation"),.fun = get_assignments, unassgined_data=unassgined_data, output=output, .parallel = parallel, testing = testing, custom = custom)  
}

