analyze <- function(modeldata, ref="actual", comp="assigned", analyzeby = "conformation"){
  require(plyr)  
  modelerror <- ddply(.dat=modeldata, .var=c(analyzeby), 
                      .fun = function(modeldata){
                        errors <- data.frame(conformation=0)
                        errors$conformation <- modeldata$conformation[1]
                        # TODO MAE, RMSE
                        errors$mae <-  mean(abs(modeldata[,comp]-modeldata[,ref]))
                        errors$rmse <- sqrt(mean((modeldata[,comp]-modeldata[,ref])^2))                        
                        errors$rmae <-  mean(abs(modeldata[,comp]-modeldata[,ref])/modeldata$error)
                        errors$rrmse <- sqrt(mean((modeldata[,comp]-modeldata[,ref])^2/(modeldata$error*modeldata$error)))
                        errors$prob <-  mean(abs(modeldata[,comp]-modeldata[,ref])<1.5*modeldata$error)                        
                        errors$r <- cor(modeldata[,comp], modeldata[,ref], method = "pearson")
                        errors$rho <- cor(modeldata[,comp], modeldata[,ref], method = "spearman")
                        errors$tau <- cor(modeldata[,comp], modeldata[,ref], method = "kendall")
                        errors
                      })
  return(modelerror)
}

analyze2 <- function(modeldata, ref="actual", comp="assigned", analyzeby = "conformation"){
  require(plyr)  
  modelerror <- ddply(.dat=modeldata, .var=analyzeby, 
                      .fun = function(modeldata){
                        errors <- data.frame(conformation=0)
                        errors$conformation <- modeldata$conformation[1]
                        # TODO MAE, RMSE
                        errors$mae <-  round(mean(abs(modeldata[,comp]-modeldata[,ref])), 2)
                        errors$rmae <-  round(mean(abs(modeldata[,comp]-modeldata[,ref])/modeldata$error), 2)
                        errors$rrmse <- sqrt(mean((modeldata[,comp]-modeldata[,ref])^2/(modeldata$error*modeldata$error)))
                        errors$r <- 1-cor(modeldata[,comp], modeldata[,ref], method = "pearson")
                        errors$rho <- 1-cor(modeldata[,comp], modeldata[,ref], method = "spearman")
                        errors$tau <- 1-cor(modeldata[,comp], modeldata[,ref], method = "kendall")                        
                        errors
                      })
  return(modelerror)
}

merge_data <- function(rna, predictor="ramsey", iters=seq(1000, 10000, 1000), conformations=1:15, weighted=FALSE, corrected=FALSE){	
	k <- 0
	for (i in conformations){
		for (iter in iters){
			suppressWarnings(try(silent=TRUE, {
			  if(weighted){
			    if(corrected){
			  		file <- paste("assignments/assigned_shifts_weighted_corrected_", predictor, "_", rna, "_", iter, "_", i, ".txt", sep ="")
			  	} else {
			  		file <- paste("assignments/assigned_shifts_weighted_uncorrected_", predictor, "_", rna, "_", iter, "_", i, ".txt", sep ="")
			  	}
			  } else {
			    if(corrected){
			  		file <- paste("assignments/assigned_shifts_unweighted_corrected_", predictor, "_", rna, "_", iter, "_", i, ".txt", sep ="")
			  	} else {
			  		file <- paste("assignments/assigned_shifts_unweighted_uncorrected_", predictor, "_", rna, "_", iter, "_", i, ".txt", sep ="")
			  	}
			  }
				tmp <- read.table(file, col.names = c("conformation", "resid", "resname", "nucleus", "actual", "assigned", "predicted"))
				tmp$iter <- iter
				if (k==0){
					cs <- tmp
					k <- 1
				} else {
					cs <- rbind(cs, tmp)
				}
			}))		
		} 	
	}
	return(cs)
}

analyze_new <- function(rna="1LDZ", predictor="ramsey", it=seq(1000, 10000, 1000), conformations=1, analyzeby=c("iter"),ref="actual", comp="assigned",  weighted=FALSE, corrected=FALSE){
	cs <- merge_data(rna, predictor, it, conformations, weighted=weighted, corrected=corrected)
	accu_file <- paste("data/", predictor,"_accuracy_resname_nucleus.txt", sep="")
	accu <- read.table(accu_file)
	colnames(accu) <- c("resname", "nucleus", "error")
	cs <- merge(cs, accu, by = c("resname", "nucleus"))
	analyze(cs, ref=ref, comp=comp, analyzeby=analyzeby)
}

summarize_rnas <- function(predictor="ramsey", start=1000, stop=10000, stride=1000, conformations=1:15, errortype="tau", rnas = unlist(strsplit("1ZC5 2LUN 2M4W 2MNC 2N6T 2M5U 2MXL 2NCI 1UUU 2JWV 2M12 2N2O 2QH2 1XHP 2K66 2LQZ 2M21 2MEQ 2N2P 2N7X 2QH4 2M22 2MFD 2N4L 1KKA 1PJY 2LI4 2LUB 4A4S 1L1W 2LBJ 2LK3 4A4T 1LC6 2FDT 2LBL 2LP9 2LV0 4A4U 1LDZ 2LDL 2LPA 1NC0 2LDT 1OW9 1YSV 2KOC 2LHP 2LU0 1SCL 1Z2J 2M24 2N6S", " "))){
	for (i in seq_along(rnas)){
		rna <- rnas[i]
		tmp <- summarize(rna, predictor, seq(start, stop, stride), conformations=conformations, errortype=errortype, ref = "actual", comp = "assigned")
		tmp$id <- rna
		if (i==1){
			summ <- tmp
		} else {
			summ <- rbind(summ, tmp)
		}			
	}
	return(summ)
}

nslr <- function(response){
	#' Sum of Logarithmic Ranks Function
	#'
	#' This function allows you to compute the normalized sum of logarithmic ranks
	#' @param response vector of 0 (inactives) and 1 (actives) that was sorted based on some scores (e.g., agreement between measured and predicted shifts)
	#' @export
	#' @examples
	#' random_nslr(sample(c(rep(0,100),rep(1,10))))
  
  ri <- which(response==1)
  N <- length(response)
  i <-  1:length(ri)
  SLRmax <- -sum(log(i/N))
  return(-sum(log(ri/N))/SLRmax)
} 

get_nslr <- function(rna="1SCL", predictor="larmord", it=seq(1,1,1), conformations=1:40, analyzeby=c("conformation"), ref="actual", comp="assigned", weighted=FALSE, corrected=FALSE){
	r <- analyze_new(rna=rna, predictor=predictor, it=seq(1,1,1), conformations=conformations, analyzeby=c("conformation"), ref=ref, comp=comp,  weighted=weighted, corrected=corrected)
	rmsd <- read.table(paste("struct_info/", rna, ".txt", sep=""), col.names=c("id", "decoy", "rmsd", "tm", "gdt", "conformation"))
	r <- merge(r, rmsd)
	r$status <- 0
	r$status[r$rmsd<2.5] <- 1
	r <- r[order(r$mae), c("conformation", "mae", "rmsd", "gdt", "status")]
  r$nslr <- nslr(r$status)
  return(r)	
}

gather_data <- function(predictor="ramsey", iteration=1, conformations=1:1, weighted=FALSE, corrected=FALSE, nuclei=c("C1'", "C2'", "C3'", "C4'", "C5'", "C2", "C5", "C6", "C8", "H1'", "H2'", "H3'", "H4'", "H5'", "H5''", "H2", "H5", "H6", "H8"), rnas = unlist(strsplit("1ZC5 2LUN 2M4W 2MNC 2N6T 2M5U 2MXL 2NCI 1UUU 2JWV 2M12 2N2O 2QH2 1XHP 2K66 2LQZ 2M21 2MEQ 2N2P 2N7X 2QH4 2M22 2MFD 2N4L 1KKA 1PJY 2LI4 2LUB 4A4S 1L1W 2LBJ 2LK3 4A4T 1LC6 2FDT 2LBL 2LP9 2LV0 4A4U 1LDZ 2LDL 2LPA 1NC0 2LDT 1OW9 1YSV 2KOC 2LHP 2LU0 1SCL 1Z2J 2M24 2N6S", " "))){
  for (i in seq_along(rnas)){
    suppressWarnings(try(silent=TRUE, {
      rna <- rnas[i]
      tmp <- merge_data(rna, predictor, iteration, conformations, weighted, corrected)
      accu_file <- paste("data/", predictor,"_accuracy_resname_nucleus.txt", sep="")
      accu <- read.table(accu_file)
      colnames(accu) <- c("resname", "nucleus", "error")
      tmp <- merge(tmp, accu, by = c("resname", "nucleus"))
      tmp$id <- rna
      tmp$diff <- tmp$assigned - tmp$actual
      tmp$wt_diff <- tmp$diff/tmp$error
      if (i==1){
        summ <- tmp
      } else {
        summ <- rbind(summ, tmp)
      }			
    }))
  }
  summ <-  summ[summ$nucleus %in% nuclei, ]
  summ <- summ[summ$actual != 0.0,]
  summ$type <- "carbon"
  summ$type[grepl("H", summ$nucleus)] <- "proton"
  summ$type[grepl("N", summ$nucleus)] <- "nitrogen"
  return(summ)
}

less_to_expected <- function(x, threshold) {mean(abs(x) < threshold)}

make_accuracy_table <- function(cs, thresholds=seq(0.50, 4.00, 0.50), round=3){
  errs <- NULL
  for (threshold in thresholds){
    errs <- c(errs, round(less_to_expected(cs$wt_diff, threshold), round))
  }
  errs <- data.frame(thresholds, errs)
  return(errs)
}

peek_at_data <- function(cs, rna, nuclei=c("C1'", "C2'", "C3'", "C4'", "C5'", "C2", "C5", "C6", "C8", "H1'", "H2'", "H3'", "H4'", "H5'", "H5''", "H2", "H5", "H6", "H8") ){
	cs <- cs[cs$id %in% rna, ]
	cs <- cs[cs$nucleus %in% nuclei, ]
	return(cs)
}

export_data_for_analysis <- function(){
	# export all the data needed for analysis
	cs_larmord_unweighted_uncorrected <- gather_data(predictor="larmord", weighted=FALSE, corrected=FALSE, conformations=1:40)
	cs_larmord_unweighted_corrected <- gather_data(predictor="larmord", weighted=FALSE, corrected=TRUE, conformations=1:40)
	cs_larmord_weighted_uncorrected <- gather_data(predictor="larmord", weighted=TRUE, corrected=FALSE, conformations=1:40)
	cs_larmord_weighted_corrected <- gather_data(predictor="larmord", weighted=TRUE, corrected=TRUE, conformations=1:40)

	cs_ramsey_unweighted_uncorrected <- gather_data(predictor="ramsey", weighted=FALSE, corrected=FALSE, conformations=1:40)
	cs_ramsey_unweighted_corrected <- gather_data(predictor="ramsey", weighted=FALSE, corrected=TRUE, conformations=1:40)
	cs_ramsey_weighted_uncorrected <- gather_data(predictor="ramsey", weighted=TRUE, corrected=FALSE, conformations=1:40)
	cs_ramsey_weighted_corrected <- gather_data(predictor="ramsey", weighted=TRUE, corrected=TRUE, conformations=1:40)
	
	data <- list(luu=cs_larmord_unweighted_uncorrected, luc=cs_larmord_unweighted_corrected, lwu=cs_larmord_weighted_uncorrected, lwc=cs_larmord_weighted_corrected, ruu=cs_ramsey_unweighted_uncorrected, ruc=cs_ramsey_unweighted_corrected, rwu=cs_ramsey_weighted_uncorrected, rwc=cs_ramsey_weighted_corrected)
	save(data, file="data/analysis_data.RData")	
}


get_rmsd_ranges <- function( rnas = unlist( strsplit("1ZC5 2LUN 2M4W 2MNC 2N6T 2M5U 2MXL 2NCI 1UUU 2JWV 2M12 2N2O 2QH2 1XHP 2K66 2LQZ 2M21 2MEQ 2N2P 2N7X 2QH4 2M22 2MFD 2N4L 1KKA 1PJY 2LI4 2LUB 4A4S 1L1W 2LBJ 2LK3 4A4T 1LC6 2FDT 2LBL 2LP9 2LV0 4A4U 1LDZ 2LDL 2LPA 1NC0 2LDT 1OW9 1YSV 2KOC 2LHP 2LU0 1SCL 1Z2J 2M24 2N6S", " ") ) ){
	for ( rna in rnas ){
		rmsd <- read.table(paste("/home/afrankz/workspace/unassigned_paper/0_assignments/struct_info/", rna, ".txt", sep=""), col.names=c("id", "decoy", "rmsd", "tm", "gdt", "conformation"))
		cat( sprintf( "%s %s %s\n", rna, nrow(rmsd), max( rmsd$rmsd ) ) )
	}
}

random_assigner <- function(cs){
	# randomly assign chemical shifts
	cs$random <- sample(cs$actual, length(cs$actual))
	return(cs)
}

randomly_assign <- function(cs){
	require(plyr)
	# apply random_assigner to a given chemical shift data set
	cs <- ddply(.dat = cs, .var = c("id", "type"), .fun =  random_assigner)
	return(cs)
}

get_accuracy <- function(cs, analyzeby = "type", ref = "assigned", comp = "actual"){
	nuclei <- c( "H1'", "H2'", "H3'", "H4'", "H5'", "H5''", "H2", "H5", "H6", "H8", "C1'", "C2'", "C3'", "C4'", "C5'", "C2", "C5", "C6", "C8")
	cs <- randomly_assign(cs)
	acc <- analyze2(cs, analyzeby = analyzeby, ref = ref, comp = comp)
	return(acc)
}


combine_accuracy_data <- function(data){  
	# combines larmord and ramsey assignment accuracy data
	tmp <- get_accuracy(cs_larmord_unweighted_uncorrected, analyzeby=c("id", "type"))
	tmp$group <- "luu"
	acc <- tmp

	tmp <- get_accuracy(cs_larmord_unweighted_corrected, analyzeby=c("id", "type"))
	tmp$group <- "luc"
	acc <- rbind(acc, tmp)

	tmp <- get_accuracy(cs_larmord_weighted_corrected, analyzeby=c("id", "type"))
	tmp$group <- "lwc"
	acc <- rbind(acc, tmp)

	tmp <- get_accuracy(cs_larmord_weighted_uncorrected, analyzeby=c("id", "type"))
	tmp$group <- "lwu"
	acc <- rbind(acc, tmp)  

	tmp <- get_accuracy(cs_ramsey_unweighted_uncorrected, analyzeby=c("id", "type"))
	tmp$group <- "ruu"
	acc <- rbind(acc, tmp)

	tmp <- get_accuracy(cs_ramsey_unweighted_corrected, analyzeby=c("id", "type"))
	tmp$group <- "ruc"
	acc <- rbind(acc, tmp)

	tmp <- get_accuracy(cs_ramsey_weighted_uncorrected, analyzeby=c("id", "type"))
	tmp$group <- "rwu"
	acc <- rbind(acc, tmp)

	tmp <- get_accuracy(cs_ramsey_weighted_corrected, analyzeby=c("id", "type"))
	tmp$group <- "rwc"
	acc <- rbind(acc, tmp)

	return(acc)
}


combine_accuracy_data <- function(data, conformer=1){  
  # combines larmord and ramsey assignment accuracy data into a flat data.frame
  # 0 - initialize object to store all results
  accu <- NULL
  # 1 - codes for each dataset
  sets <- c("luu", "luc", "lwu", "lwc", "ruu", "ruc", "rwu", "rwc")
  # 2 - loop over data set
  for (set in sets){
    # 2-0 - get accuracy stats
    data_tmp <- data[[set]]
    data_tmp <- subset(data_tmp, conformation==conformer)
    tmp <- get_accuracy(data_tmp, analyzeby=c("id", "type"))
    # 2-1 - set group id to dataset code
    tmp$group <- set
    # 2-2 - combine accuracy data
    ifelse(is.null(accu), accu <- tmp, accu <- rbind(accu, tmp))
  }
  return(accu)
}


get_boxplot <- function(acc, metric = "rmae", ty = "proton", cols = c(rep("black"), 4), codes = c("luu", "lwu", "ruu", "rwu"), lim =c(0, 1), names = NULL, ylab = "", plot=TRUE){
  # make boxplots
  acc <- subset(acc [acc$group %in% codes, ], type == ty)
  acc$metric <- acc[, metric]
  p <- boxplot(metric~type+group, acc, ylim = lim, las=2, col = cols, names = names, ylab = ylab, plot = plot)
  return(p)
}


get_accuracy_subset <- function(acc, grp = "luu", ty = "proton", sortby = "mae"){
  # return subset of the accuracy info sorted by specified column
  acc <- subset(acc, type==ty & group==grp)
  acc <- acc[order(acc[, sortby]), ]
  return(acc)
}

get_data_subset <- function(data, rna, grp = "luu"){
  # return subset of the data based on rna and dataset code
  data <- data[[grp]]
  data <- subset(data, id==rna)
  return(data)
}


pyshifts_ready_data <- function(data, rna, set="luu", filename="none.txt", cs="predicted", conformer=1, predicted_file = TRUE){
  tmp <- data[[set]]
  tmp <- subset(tmp, conformation==conformer&id==rna)
  if (predicted_file){
    tmp <- tmp[, c("conformation", "resid", "resname", "nucleus", cs, "id")]
  } else {
    tmp <- tmp[, c("resname", "resid", "nucleus", cs, "error")]
  }
  write.table(tmp, file = filename, col.names = FALSE, row.names = FALSE, quote = FALSE)
}


evaluate_nslr <- function(data, rna, ref="actual", comp="assigned", analyzeby=c("conformation"), sortby = "mae", struct_location="/Users/atfrank/Desktop/summer-papers/analysis/struct_info"){
  # returns NSLR along with errors for each conformation
  r <- analyze2(data, ref=ref, comp=comp, analyzeby=analyzeby)
  rmsd <- read.table(paste(struct_location, "/", rna, ".txt", sep=""), col.names=c("id", "decoy", "rmsd", "tm", "gdt", "conformation"))
  r <- merge(r, rmsd)
  r$status <- 0
  r$status[r$rmsd<2.5] <- 1
  r <- r[order(r[, sortby]), c("conformation", sortby, "rmsd", "gdt", "status")]
  nslr <- nslr(r$status)
  return(list(errors=r, nslr=nslr))
}

collect_nslrs <- function(data, rnas = unlist(strsplit("1ZC5 2LUN 2M4W 2MNC 2N6T 2M5U 2MXL 2NCI 1UUU 2JWV 2M12 2N2O 2QH2 1XHP 2K66 2LQZ 2M21 2MEQ 2N2P 2N7X 2QH4 2M22 2MFD 2N4L 1KKA 1PJY 2LI4 2LUB 4A4S 1L1W 2LBJ 2LK3 4A4T 1LC6 2FDT 2LBL 2LP9 2LV0 4A4U 1LDZ 2LDL 2LPA 1NC0 2LDT 1OW9 1YSV 2KOC 2LHP 2LU0 1SCL 1Z2J 2M24 2N6S", " ")), grp = "luu", ref = "actual", comp = "assigned", analyzeby = c("conformation"), sortby = "rmae"){
  # collects the NSLRs for all RNAs using a given error metric and for a specific dataset specified by grp code
  nslr <- NULL
  for (rna in rnas){
    tmp <- get_data_subset(data, rna, grp = grp)
    tmp <- evaluate_nslr(tmp, rna = rna, ref = ref, comp = comp, analyzeby = analyzeby, sortby = sortby)
    tmp <- tmp$nslr
    ifelse(is.null(nslr), nslr <- tmp, nslr <- c(nslr, tmp))
  }
  return(nslr)
}

in_training_sets <- function(rna, larmord_train_rnas = unlist(strsplit("1KKA 1L1W 1LC6 1LDZ 1NC0 1OW9 1PJY 1R7W 1R7Z 1YSV 2FDT 2GM0 2JXQ 2JXS 2K3Z 2K41 2KOC 2KYD 2LBJ 2LBL 2LDL 2LDT 2LHP 2LI4 2LK3 2LP9 2LPA 2LU0 2LUB 2LV0 2RN1 2Y95 4A4S 4A4T 4A4U", " ")), ramsey_train_rnas =  unlist(strsplit("1LC6 1LDZ 1NC0 1OW9 1PJY 1R7W 1R7Z 1UUU 1YSV 1Z2J 2GM0 2JTP 2JXQ 2JXS 2K3Z 2K41 2KOC 2L3E 2LDL", " "))){
  rna %in% unique(c(larmord_train_rnas, ramsey_train_rnas))
}


get_nslr_matrix <- function(data,  grp = "luu", ref = "actual", comp = "assigned", errortypes = c("rmae","rrmse", "r", "rho", "tau")){
  # collect nlsr matrix for ref vs. comp for a given dataset
  rnas <-  unlist(strsplit("1ZC5 2LUN 2M4W 2MNC 2N6T 2M5U 2MXL 2NCI 1UUU 2JWV 2M12 2N2O 2QH2 1XHP 2K66 2LQZ 2M21 2MEQ 2N2P 2QH4 2M22 2MFD 2N4L 1KKA 1PJY 2LI4 2LUB 4A4S 1L1W 2LBJ 2LK3 4A4T 1LC6 2FDT 2LBL 2LP9 2LV0 4A4U 1LDZ 2LDL 2LPA 1NC0 2LDT 1OW9 1YSV 2KOC 2LHP 2LU0 1SCL 1Z2J 2M24 2N6S", " "))
  nslr_matrix <- NULL
  for (errortype in errortypes){
    nslr <- collect_nslrs(data=data, rnas = rnas, grp = grp, ref = ref, comp = comp, analyzeby = c("conformation"), sortby = errortype)
    ifelse(is.null(nslr_matrix), nslr_matrix <- nslr, nslr_matrix <- cbind(nslr_matrix, nslr))
  }
  colnames(nslr_matrix) <- errortypes
  rownames(nslr_matrix) <- rnas
  return(nslr_matrix)
}


collect_errors <- function(data, grp = "luu", ref = "actual", comp = "assigned", analyzeby = c("conformation"), sortby = "rmae", rnas = unlist(strsplit("1ZC5 2LUN 2M4W 2MNC 2N6T 2M5U 2MXL 2NCI 1UUU 2JWV 2M12 2N2O 2QH2 1XHP 2K66 2LQZ 2M21 2MEQ 2N2P 2QH4 2M22 2MFD 2N4L 1KKA 1PJY 2LI4 2LUB 4A4S 1L1W 2LBJ 2LK3 4A4T 1LC6 2FDT 2LBL 2LP9 2LV0 4A4U 1LDZ 2LDL 2LPA 1NC0 2LDT 1OW9 1YSV 2KOC 2LHP 2LU0 1SCL 1Z2J 2M24 2N6S", " "))){
  # collects the NSLRs for all RNAs using a given error metric and for a specific dataset specified by grp code
  errors <- NULL
  for (rna in rnas){
    tmp <- get_data_subset(data, rna, grp = grp)
    tmp <- evaluate_nslr(tmp, rna = rna, ref = ref, comp = comp, analyzeby = analyzeby, sortby = sortby)
    tmp <- tmp$errors
    tmp$id <- rna
    ifelse(in_training_sets(rna), tmp$set <- "training", tmp$set <- "testing" )
    ifelse(is.null(errors), errors <- tmp, errors <- rbind(errors, tmp))
  }
  return(errors)
}

get_rmsd_matrix <- function(data, N=1, grps = c("luu", "luc", "ruu", "ruc"), ref = "actual", comp = "assigned", errortype = c("rmae"), rnas = unlist(strsplit("1ZC5 2LUN 2M4W 2MNC 2N6T 2M5U 2MXL 2NCI 1UUU 2JWV 2M12 2N2O 2QH2 1XHP 2K66 2LQZ 2M21 2MEQ 2N2P 2QH4 2M22 2MFD 2N4L 1KKA 1PJY 2LI4 2LUB 4A4S 1L1W 2LBJ 2LK3 4A4T 1LC6 2FDT 2LBL 2LP9 2LV0 4A4U 1LDZ 2LDL 2LPA 1NC0 2LDT 1OW9 1YSV 2KOC 2LHP 2LU0 1SCL 1Z2J 2M24 2N6S", " "))){
  # collect RMSD matrix 
  rmsd_matrix <- NULL
  for (grp in grps){
    rmsd <- collect_errors(data=data, grp = grp, ref = ref, comp = comp, analyzeby = c("conformation"), sortby = errortype, rnas = rnas)
    rmsd <- plyr::ddply(.data = rmsd, .variables = c("id"), .fun = function(x) {mean(head(x$rmsd, N))})
    rownames(rmsd) <- rmsd$id
    names <- rownames(rmsd)
    rmsd <- rmsd$V1
    ifelse(is.null(rmsd_matrix), rmsd_matrix <- rmsd, rmsd_matrix <- cbind(rmsd_matrix, rmsd))
  }
  
  colnames(rmsd_matrix) <- grps
  rownames(rmsd_matrix) <- names
  return(rmsd_matrix)
}


get_rna_sets <- function(rnas = unlist(strsplit("1ZC5 2LUN 2M4W 2MNC 2N6T 2M5U 2MXL 2NCI 1UUU 2JWV 2M12 2N2O 2QH2 1XHP 2K66 2LQZ 2M21 2MEQ 2N2P 2QH4 2M22 2MFD 2N4L 1KKA 1PJY 2LI4 2LUB 4A4S 1L1W 2LBJ 2LK3 4A4T 1LC6 2FDT 2LBL 2LP9 2LV0 4A4U 1LDZ 2LDL 2LPA 1NC0 2LDT 1OW9 1YSV 2KOC 2LHP 2LU0 1SCL 1Z2J 2M24 2N6S", " ")), larmord_train_rnas = unlist(strsplit("1KKA 1L1W 1LC6 1LDZ 1NC0 1OW9 1PJY 1R7W 1R7Z 1YSV 2FDT 2GM0 2JXQ 2JXS 2K3Z 2K41 2KOC 2KYD 2LBJ 2LBL 2LDL 2LDT 2LHP 2LI4 2LK3 2LP9 2LPA 2LU0 2LUB 2LV0 2RN1 2Y95 4A4S 4A4T 4A4U", " ")), ramsey_train_rnas =  unlist(strsplit("1LC6 1LDZ 1NC0 1OW9 1PJY 1R7W 1R7Z 1UUU 1YSV 1Z2J 2GM0 2JTP 2JXQ 2JXS 2K3Z 2K41 2KOC 2L3E 2LDL", " "))){
  ramsey <- rnas[(rnas %in% ramsey_train_rnas)]
  rest <- rnas[!(rnas %in% ramsey_train_rnas)]
  larmord <- rest[(rest %in% larmord_train_rnas)]
  rest <- rest[!(rest %in% larmord_train_rnas)]
  larmord_train_rnas <- rnas[rnas %in% larmord_train_rnas]
  ramsey_train_rnas <- rnas[rnas %in% ramsey_train_rnas]
  
  return(list(all=c(rest, larmord, ramsey), ramsey=ramsey, larmord=larmord, larmord_train_rnas = larmord_train_rnas, ramsey_train_rnas = ramsey_train_rnas))
}

make_levelplot <- function(rmsd_matrix, at = seq(0, 7, 0.5), cols = c("red", "blue")){
  # make levelplot (heatmap) from rmsd_matrices
  require(lattice)
  
  myPanel <- function(x, y, z, ...) {
    panel.levelplot(x,y,z,...)
    panel.grid(h = 51, v = 3, col = "black")
    # panel.text(x, y, round(z,1))
  }
  
  rnas <- get_rna_sets()
  rmsd_matrix <- rmsd_matrix[rnas$all, ]
  colRamp <- colorRampPalette(cols[1:2])
  grid <- expand.grid(X=1:ncol(rmsd_matrix), Y=rownames(rmsd_matrix))
  grid$Z <- as.vector(t(rmsd_matrix))
  levelplot(Z ~ X*Y, grid, panel = myPanel, at=at, outer = FALSE, las = 2, pretty= TRUE, useRaster = TRUE, colorkey = TRUE, col.regions = colRamp, xlab = "", ylab = "")
}

save_levelplots <- function(matrix, fname = "test.pdf", width = 2.30, height = 7.27, at=seq(0, 7, 1), cols = c("white","black")){
  pdf(file=fname, width = width, height = height)
  make_levelplot(matrix, cols = cols, at = at)
  dev.off()
}


fraction_lower_than <- function(data, at = seq(1, 7, 0.5)){
  # get fraction of test tests with RMSD than given thresholds
  fraction <- NULL
  for (a in at){
    tmp <- data < a
    tmp <- colMeans(tmp)
    ifelse(is.null(fraction), fraction <- tmp, fraction <- rbind(fraction, tmp))
  }
  rownames(fraction) <- at
  return(fraction)
}

plot_fraction_lower_than <- function(data, set="all", larmord_train_rnas = unlist(strsplit("1KKA 1L1W 1LC6 1LDZ 1NC0 1OW9 1PJY 1R7W 1R7Z 1YSV 2FDT 2GM0 2JXQ 2JXS 2K3Z 2K41 2KOC 2KYD 2LBJ 2LBL 2LDL 2LDT 2LHP 2LI4 2LK3 2LP9 2LPA 2LU0 2LUB 2LV0 2RN1 2Y95 4A4S 4A4T 4A4U", " ")), ramsey_train_rnas =  unlist(strsplit("1LC6 1LDZ 1NC0 1OW9 1PJY 1R7W 1R7Z 1UUU 1YSV 1Z2J 2GM0 2JTP 2JXQ 2JXS 2K3Z 2K41 2KOC 2L3E 2LDL", " "))){
  # function that plots fraction of dataset with RMSD lower than thresholds
  # separate larmord and ramsey data
  larmord <- data[,c("luu", "luc")]
  ramsey <- data[,c("ruu", "ruc")]
  
  # select subset of data corresponding to the dataset requested
  rnas_sets <- get_rna_sets()
  larmord_train_rnas <- rnas_sets$all[rnas_sets$all %in%larmord_train_rnas]
  ramsey_train_rnas <- rnas_sets$all[rnas_sets$all %in%ramsey_train_rnas]
  
  # if training set
  if(set == "training"){
    larmord <- larmord[larmord_train_rnas, ]
    ramsey <- ramsey[ramsey_train_rnas, ]
  }
  # if testing set
  if(set == "testing"){
    rnas <- rnas_sets$all[!(rnas_sets$all %in% larmord_train_rnas)]
    larmord <- larmord[rnas, ]
    rnas <- rnas_sets$all[!(rnas_sets$all %in% ramsey_train_rnas)]
    ramsey <- ramsey[rnas, ]
  }
  
  rmsd_larmord <- colMeans(larmord)
  rmsd_ramsey <- colMeans(ramsey)
  
  larmord <- fraction_lower_than(larmord)
  ramsey <- fraction_lower_than(ramsey)
  
  xs <- as.numeric(rownames(larmord))
  
  plot(x = xs, y = larmord[,1], type = "l", lwd = 2, col = "black", ylim = c(0, 1), ylab = "", xlab = "")
  abline(v=mean(rmsd_larmord[1]), lwd="2", col = "black")
  lines(x = xs, y = larmord[,2], type = "l", lwd = 2, lty=2, col = "black")
  abline(v=mean(rmsd_larmord[2]), lwd="2", col = "black", lty=2)
  lines(x = xs, y = ramsey[,1], type = "l", lwd = 2, col = "gray")
  abline(v=mean(rmsd_ramsey[1]), lwd="2", col = "gray")
  lines(x = xs, y = ramsey[,2], type = "l", lwd = 2, lty=2, col = "gray")
  abline(v=mean(rmsd_ramsey[2]), lwd="2", col = "gray", lty=2)
  abline(v=2.5, lwd="3", col = "red")
  return(list(larmord_fraction = larmord,  larmord_rmsd=rmsd_larmord, ramsey_rmsd=rmsd_ramsey, ramsey_fraction = ramsey))
}


get_median_and_iqr <- function(p_p_1){
  p <- p_p_1$stats
  iqr <- p[4,]-p[2,]
  median <- p[3,] 
  return(list(median=median, iqr=iqr))
}


summarize_accuracy_for_paper <- function(fname="proton_accuracy.pdf", nuclei = "proton", codes = c("luu", "ruu"),  cols = c("white", "grey"), lim = c(0, 1), width = 2.36, height = 5.61){
  rnas <- get_rna_sets()
  acc$in_larmord_training_set <- FALSE
  acc$in_larmord_training_set[acc$id %in% rnas$larmord_train_rnas] <- TRUE
  
  acc$in_ramsey_training_set <- FALSE
  acc$in_ramsey_training_set[acc$id %in% rnas$ramsey_train_rnas] <- TRUE
  pdf(file=fname, width = 2.36, height = 5.61)
  p_p_1 <- get_boxplot(acc, metric = "mae", ty = nuclei, codes = codes, lim = lim, cols = cols, names = c("LARMORD", "RAMSEY"), ylab = "assignment error (ppm)")
  print(p_p_1)
  dev.off()
  p_p_1_l_train <- get_boxplot(subset(acc, in_larmord_training_set == TRUE), plot=FALSE, metric = "mae", ty = nuclei, codes = codes[1], lim = c(0, 1), cols = c("white"), names = c("LARMORD"), ylab = "assignment error (ppm)")
  p_p_1_l_test <- get_boxplot(subset(acc, in_larmord_training_set == FALSE), plot=FALSE, metric = "mae", ty = nuclei, codes = codes[1], lim = c(0, 1), cols = c("white"), names = c("LARMORD"), ylab = "assignment error (ppm)")
  p_p_1_r_train <- get_boxplot(subset(acc, in_ramsey_training_set == TRUE),  plot=FALSE, metric = "mae", ty = nuclei, codes = codes[2], lim = c(0, 1), cols = c("white"), names = c("RAMSEY"), ylab = "assignment error (ppm)")
  p_p_1_r_test <- get_boxplot(subset(acc, in_ramsey_training_set == FALSE),  plot=FALSE, metric = "mae", ty = nuclei, codes = codes[2], lim = c(0, 1), cols = c("white"), names = c("RAMSEY"), ylab = "assignment error (ppm)")
  
  r <- get_median_and_iqr(p_p_1)
  cat(sprintf("For all %ss, LARMORD/RAMSEY median error:      %4.2f/%4.2f\n", nuclei, round(r$median[1], 2), round(r$median[2], 2)))
  cat(sprintf("For all %ss, LARMORD/RAMSEY IQR:               %4.2f/%4.2f\n", nuclei, round(r$iqr[1], 2), round(r$iqr[2], 2)))
  cat(sprintf("\n"))

  ltrain <- get_median_and_iqr(p_p_1_l_train)
  rtrain <- get_median_and_iqr(p_p_1_r_train)
  cat(sprintf("For training %ss, LARMORD/RAMSEY median error: %4.2f/%4.2f\n", nuclei, round(ltrain$median[1], 2), round(rtrain$median[1], 2)))
  cat(sprintf("For training %ss, LARMORD/RAMSEY IQR:          %4.2f/%4.2f\n", nuclei, round(ltrain$iqr[1], 2), round(rtrain$iqr[1], 2)))
  cat(sprintf("\n"))
  
  ltest <- get_median_and_iqr(p_p_1_l_test)
  rtest <- get_median_and_iqr(p_p_1_r_test)
  cat(sprintf("For testing %ss, LARMORD/RAMSEY median error:  %4.2f/%4.2f\n", nuclei, round(ltest$median[1], 2), round(rtest$median[1], 2)))
  cat(sprintf("For testing %ss, LARMORD/RAMSEY IQR:           %4.2f/%4.2f\n", nuclei, round(ltest$iqr[1], 2), round(rtest$iqr[1], 2)))
  cat(sprintf("\n"))
}

rnawise_shifts <- function(tmp, nuclei=c("C1'","C2'","C3'","C4'","C5'","C2","C5","C6","C8")){
  # helper function used to output residue-wise chemical shifts
  shifts <- NULL
  for (nucleus in nuclei){
    s <- abs(round(tmp$cs[tmp$nucleus==nucleus], 2))
    shifts <- c(shifts,ifelse(!is.null(s), s, NA))
  }
  return(shifts)
}

load_detected_referencing_errors <- function(filename){
  # loads detected referencing error data
  ref_errors <- read.table(filename, col.names = c("method", "id", "cs",  "nucleus"))
  r <- get_rna_sets()
  return(ref_errors[ref_errors$id %in% r$all,])
}

make_ref_error_table <- function(filename = "detected_errors.txt", nuclei = c("C1'","C2'","C3'","C4'","C5'","C2","C5","C6","C8")){
  # makes a table of referencing errors suitable for manuscript
  ref_errors <- load_detected_referencing_errors(filename)
  ref_errors <- plyr::ddply(.data = ref_errors, .variables = c("method", "id"), .fun = rnawise_shifts, nuclei)
  colnames(ref_errors) <- c("method", "id", nuclei)
  return(ref_errors)
}

make_timing_plot <- function(infname = "timinings.txt", outfname = "timing.pdf", width = 4.41, height = 3.96){
  # plot of the time required for assignment per conformer
  pdf(file=outfname, width = width, height = height)
  timing <- read.table(infname, col.names = c("id", "n_actual", "n_pred", "conformers", "start", "end"))
  bmrb <- unique(read.table("pdb_bmrb.txt", header=T))
  r <- get_rna_sets()
  timing <- merge(timing, bmrb, all = TRUE, by = c("id"))
  timing <- timing[timing$id %in% r$all, ]
  rownames(timing) <- timing$id
  timing$duration <- as.difftime(as.character(timing$end), units = "secs")-as.difftime(as.character(timing$start), units = "secs")
  timing$duration2 <- as.double(as.difftime(as.character(timing$end), units = "secs")-as.difftime(as.character(timing$start), units = "secs"))
  timing$per_conformer <- timing$duration2/timing$conformers
  timing <- timing[order(timing$n_pred), ]
  plot(timing$n_pred, timing$per_conformer, type = "h", lwd = "3")
  abline(h=mean(timing$per_conformer), col = "black", lwd = "2")
  abline(h=mean(timing$per_conformer), col = "red", lwd = "2", lty = 2)
  dev.off()
  cat(sprintf("mean assignment time per conformer: %4.3f\n", mean(timing$per_conformer)))
  timing <- timing[r$all, ]
  timing$index <- 1:nrow(timing)
  timing <- timing[order(timing$index, decreasing = TRUE), ]
  return(timing)
}

perfect_assignments <- function(data, conformer=1, analyzeby=c("id")){
  # check how many perfect assignments were made for each RNA
  require(plyr)
  data_tmp <- subset(data, conformation==conformer)
  plyr::ddply(.data = data_tmp, .variables = analyzeby, .fun = count_perfect_assignments)
}

count_perfect_assignments <- function(data){
  # count how many perfect assignments for a given
  n_perfect <- sum(data$diff==0)
  fraction_perfect <- mean(data$diff==0)
  n_less_than_error <- sum(abs(data$diff) < abs(data$error))
  fraction_less_than_error <- mean(abs(data$diff) < abs(data$error))
  return(data.frame(n_perfect=n_perfect, fraction_perfect=fraction_perfect, n_less_than=n_less_than_error, fraction_less_than=fraction_less_than_error))
}

summarize_perfect_assignment_accuracy <- function(data, sets=c("all", "training", "testing"), grps = c("luc", "ruc"), larmord_train_rnas = unlist(strsplit("1KKA 1L1W 1LC6 1LDZ 1NC0 1OW9 1PJY 1R7W 1R7Z 1YSV 2FDT 2GM0 2JXQ 2JXS 2K3Z 2K41 2KOC 2KYD 2LBJ 2LBL 2LDL 2LDT 2LHP 2LI4 2LK3 2LP9 2LPA 2LU0 2LUB 2LV0 2RN1 2Y95 4A4S 4A4T 4A4U", " ")), ramsey_train_rnas =  unlist(strsplit("1LC6 1LDZ 1NC0 1OW9 1PJY 1R7W 1R7Z 1UUU 1YSV 1Z2J 2GM0 2JTP 2JXQ 2JXS 2K3Z 2K41 2KOC 2L3E 2LDL", " "))){
  # returns a summarize of the analysis of fraction of assignments that were perfect and within the accuracy of the predictors
  rnas_sets <- get_rna_sets()
  larmord_train_rnas <- rnas_sets$all[rnas_sets$all %in%larmord_train_rnas]
  ramsey_train_rnas <- rnas_sets$all[rnas_sets$all %in%ramsey_train_rnas]
  
  # initialize object that will store final results
  fractions <- NULL
  
  for (i in seq_along(sets)){
    set <- sets[i]
    for (j  in seq_along(grps)){
      grp <- grps[j]
      rnas <- rnas_sets$all
      if(set == "training"){
        if (grp %in% c("luu","luc")){
          rnas <- rnas_sets$all[(rnas_sets$all %in% larmord_train_rnas)]
        } else {
          rnas <- rnas_sets$all[(rnas_sets$all %in% ramsey_train_rnas)]
        }
      }
      # if testing set
      if(set == "testing"){
        if (grp %in% c("luu","luc")){
          rnas <- rnas_sets$all[!(rnas_sets$all %in% larmord_train_rnas)]
        } else {
          rnas <- rnas_sets$all[!(rnas_sets$all %in% ramsey_train_rnas)]
        }
      }
      
      # get fraction 
      fraction <- perfect_assignments(data[[grp]])
      fraction$set <- set
      fraction$grp <- grp
      fraction <- fraction[(fraction$id %in% rnas), ]
      ifelse(is.null(fractions), fractions <- fraction, fractions <- rbind(fractions, fraction))
      cat(sprintf("%s %s %4.3f %4.3f\n", set, grp, median(fraction$fraction_perfect), median(fraction$fraction_less_than)))
    }
  }
  return(fractions)
}

summarize_identification_accuracy <- function(native_rmsd_cutoff = "2.5", matrices = c("actual_v_assigned_rmsd_1", "assigned_v_predicted_rmsd_1", "actual_v_predicted_rmsd_1")){
  # function that produces all the plots needed for summarizing the ability to identify native-like structures
  width <- 4.6
  height <- 4.0
  for (i in seq_along(matrices)){
    for (type in c("all", "training", "testing")){
      data <- get(matrices[i])
      if (!is.null(data)){
        pdf(file=paste("fraction_lower_than_", type, "_", matrices[i], ".pdf", sep = ""), width = width, height = height)
        results <- plot_fraction_lower_than(data, set = type)
        dev.off()
        cat(sprintf(" Mean RMSD for LARMORD and RAMSEY %s %s: %4.2f/%4.2f & %4.2f/%4.2f\n", type, matrices[i], mean(results$larmord_rmsd[1]), mean(results$larmord_rmsd[2]), mean(results$ramsey_rmsd[1]), mean(results$ramsey_rmsd[2])))
        cat(sprintf(" Fraction less than %s %s %s           : %4.2f/%4.2f & %4.2f/%4.2f\n", native_rmsd_cutoff, type, matrices[i], mean(results$larmord_fraction[native_rmsd_cutoff, 1]), mean(results$larmord_fraction[native_rmsd_cutoff, 2]), mean(results$ramsey_fraction[native_rmsd_cutoff, 1]), mean(results$ramsey_fraction[native_rmsd_cutoff, 2])))
        
      } else {
        print(nrow(data))
        cat(sprintf("%s object not found\n", matrices[i]))
      }
    }
    cat(sprintf("\n"))
  }
}


calculate_errors_structure_correlation <- function(tmp, structural_similarity_measure = "rmsd", errortype="rmae"){
  # calculates the correlation between the RMAE and RMSD
  R = cor(tmp[, structural_similarity_measure], tmp[, errortype], method="pear")
  rho = cor(tmp[, structural_similarity_measure], tmp[, errortype], method="spear")
  tau = cor(tmp[, structural_similarity_measure], tmp[, errortype], method="kendall")
  return(data.frame(R=R, rho=rho, tau=tau))
}


plot_rmsd_correlations <- function(data, corel_type =  "R", col = "black", set="all", grp = "luu", addplot = FALSE, larmord_train_rnas = unlist(strsplit("1KKA 1L1W 1LC6 1LDZ 1NC0 1OW9 1PJY 1R7W 1R7Z 1YSV 2FDT 2GM0 2JXQ 2JXS 2K3Z 2K41 2KOC 2KYD 2LBJ 2LBL 2LDL 2LDT 2LHP 2LI4 2LK3 2LP9 2LPA 2LU0 2LUB 2LV0 2RN1 2Y95 4A4S 4A4T 4A4U", " ")), ramsey_train_rnas =  unlist(strsplit("1LC6 1LDZ 1NC0 1OW9 1PJY 1R7W 1R7Z 1UUU 1YSV 1Z2J 2GM0 2JTP 2JXQ 2JXS 2K3Z 2K41 2KOC 2L3E 2LDL", " "))){
  # plot the correlation between RMSD and some error metric
  rnas_sets <- get_rna_sets()
  larmord_train_rnas <- rnas_sets$all[rnas_sets$all %in%larmord_train_rnas]
  ramsey_train_rnas <- rnas_sets$all[rnas_sets$all %in%ramsey_train_rnas]
  
  # get correct rna ids
  rnas <- rnas_sets$all
  if(set == "training"){
    if (grp %in% c("luu","luc")){
      rnas <- rnas_sets$all[(rnas_sets$all %in% larmord_train_rnas)]
    } else {
      rnas <- rnas_sets$all[(rnas_sets$all %in% ramsey_train_rnas)]
    }
  }
  # if testing set
  if(set == "testing"){
    if (grp %in% c("luu","luc")){
      rnas <- rnas_sets$all[!(rnas_sets$all %in% larmord_train_rnas)]
    } else {
      rnas <- rnas_sets$all[!(rnas_sets$all %in% ramsey_train_rnas)]
    }
  }
  
  # calculate errors
  errors <- collect_errors(data, grp = grp)
  # get correlations
  correlations <- plyr::ddply(.data = errors, .variables = c("id"), .fun = calculate_errors_structure_correlation)
  rownames(correlations) <- correlations$id 
  correlations <- correlations[rnas, ]
  corel <- density(correlations[, corel_type], from = -1, to = 1)
  corelnorm <- corel$y/sum(corel$y)
  # plot correlations
  if (!addplot){
    plot(x = corel$x, y = corelnorm,  col = col, lwd = "2", type = "l")
  } else {
    lines(x = corel$x, y = corelnorm, col = col, lwd = "2")
  }
  abline(v=median(correlations[, corel_type]), lwd = "2", col = col, lty = 2)
  return(rnas)
}

collect_correlations <- function(data, set = "all", ref = "actual", comp = "comp", grps = c("luu", "luc", "ruu", "ruc"), larmord_train_rnas = unlist(strsplit("1KKA 1L1W 1LC6 1LDZ 1NC0 1OW9 1PJY 1R7W 1R7Z 1YSV 2FDT 2GM0 2JXQ 2JXS 2K3Z 2K41 2KOC 2KYD 2LBJ 2LBL 2LDL 2LDT 2LHP 2LI4 2LK3 2LP9 2LPA 2LU0 2LUB 2LV0 2RN1 2Y95 4A4S 4A4T 4A4U", " ")), ramsey_train_rnas =  unlist(strsplit("1LC6 1LDZ 1NC0 1OW9 1PJY 1R7W 1R7Z 1UUU 1YSV 1Z2J 2GM0 2JTP 2JXQ 2JXS 2K3Z 2K41 2KOC 2L3E 2LDL", " "))){
  # collect correlations for various analysis groups

  # plot the correlation between RMSD and some error metric
  rnas_sets <- get_rna_sets()
  larmord_train_rnas <- rnas_sets$all[rnas_sets$all %in%larmord_train_rnas]
  ramsey_train_rnas <- rnas_sets$all[rnas_sets$all %in%ramsey_train_rnas]
    
  correlations <- NULL
  for (grp in grps){
    # get correct rna ids
    rnas <- rnas_sets$all
    if(set == "training"){
      if (grp %in% c("luu","luc")){
        rnas <- rnas_sets$all[(rnas_sets$all %in% larmord_train_rnas)]
      } else {
        rnas <- rnas_sets$all[(rnas_sets$all %in% ramsey_train_rnas)]
      }
    }
    # if testing set
    if(set == "testing"){
      if (grp %in% c("luu","luc")){
        rnas <- rnas_sets$all[!(rnas_sets$all %in% larmord_train_rnas)]
      } else {
        rnas <- rnas_sets$all[!(rnas_sets$all %in% ramsey_train_rnas)]
      }
    }
    
    errors <- collect_errors(data, grp = grp, ref = ref, comp = comp)
    correlation <- plyr::ddply(.data = errors, .variables = c("id"), .fun = calculate_errors_structure_correlation)
    rownames(correlation) <- correlation$id 
    correlation <- correlation[rnas, ]
    
    correlation$grp <- grp
    ifelse(is.null(correlations), correlations <- correlation, correlations <- rbind(correlations, correlation))
  }
  p <- boxplot(R~grp, data=correlations, col = c("white", "grey"), ylim = c(-1, 1))
  print(p)
  return(p)
}


save_correlation_boxplots <- function(data, ref = "actual", comp = "comp"){
  width <- 2.6
  height <- 3.4
  pdf(file="correlations_all.pdf", width = width, height = height)
  correlations <- collect_correlations(data, ref = ref, comp = comp, set = "all", grps = c("luc", "ruc"))
  dev.off()
  
  pdf(file="correlations_training.pdf", width = width, height = height)
  correlations <- collect_correlations(data, ref = ref, comp = comp, set = "training", grps = c("luc", "ruc"))
  dev.off()
  
  pdf(file="correlations_testing.pdf", width = width, height = height)
  correlations <- collect_correlations(data, ref = ref, comp = comp, set = "testing", grps = c("luc", "ruc"))
  dev.off()
}

get_sensitivity <- function(tmp, rosetta_scores = FALSE){
  require(pROC)
  require(plyr)
  # return the auc
  tmp <- tmp[order(tmp$rmae), ]
  tmp$status <- as.factor(tmp$status)
  proc <- roc(response = tmp$status, predictor = tmp$rmae, plot = FALSE, c)
  nslr <- nslr(response = tmp$status)
  if(rosetta_scores){
    tmp <- tmp[order(tmp$score), ]
    tmp$status <- as.factor(tmp$status)
    proc_rosetta <- roc(response = tmp$status, predictor = tmp$score, plot = FALSE, c)
    nslr_rosetta <- nslr(response = tmp$status)
    return(data.frame(roc = proc$auc, nslr = nslr, rosetta_roc = proc_rosetta$auc, nslr_rosetta = nslr_rosetta))
  } else {
    return(data.frame(roc = proc$auc, nslr = nslr))
  }
  
}

get_sensivities <- function(data, ref = "assigned", comp = "predicted", sets=c("all", "training", "testing"), grps = c("luc", "ruc"), rosetta_scores = FALSE, larmord_train_rnas = unlist(strsplit("1KKA 1L1W 1LC6 1LDZ 1NC0 1OW9 1PJY 1R7W 1R7Z 1YSV 2FDT 2GM0 2JXQ 2JXS 2K3Z 2K41 2KOC 2KYD 2LBJ 2LBL 2LDL 2LDT 2LHP 2LI4 2LK3 2LP9 2LPA 2LU0 2LUB 2LV0 2RN1 2Y95 4A4S 4A4T 4A4U", " ")), ramsey_train_rnas =  unlist(strsplit("1LC6 1LDZ 1NC0 1OW9 1PJY 1R7W 1R7Z 1UUU 1YSV 1Z2J 2GM0 2JTP 2JXQ 2JXS 2K3Z 2K41 2KOC 2L3E 2LDL", " "))){
  # calculate the AUC from the ROC and NSLR
  rnas_sets <- get_rna_sets()
  larmord_train_rnas <- rnas_sets$all[rnas_sets$all %in%larmord_train_rnas]
  ramsey_train_rnas <- rnas_sets$all[rnas_sets$all %in%ramsey_train_rnas]
  
  # initialize object that will store final results
  sensivities <- NULL
  
  # get rosetta scores
  if(rosetta_scores){scores <- merge_rosetta_scores()}
  
  for (i in seq_along(sets)){
    set <- sets[i]
    for (j  in seq_along(grps)){
      grp <- grps[j]
      rnas <- rnas_sets$all
      if(set == "training"){
        if (grp %in% c("luu","luc")){
          rnas <- rnas_sets$all[(rnas_sets$all %in% larmord_train_rnas)]
        } else {
          rnas <- rnas_sets$all[(rnas_sets$all %in% ramsey_train_rnas)]
        }
      }
      # if testing set
      if(set == "testing"){
        if (grp %in% c("luu","luc")){
          rnas <- rnas_sets$all[!(rnas_sets$all %in% larmord_train_rnas)]
        } else {
          rnas <- rnas_sets$all[!(rnas_sets$all %in% ramsey_train_rnas)]
        }
      }
      
      # calculate errors
      errors <- collect_errors(data, grp = grp, ref = ref, comp = comp)
      if(rosetta_scores){errors <- merge(errors, scores, by = c("id", "conformation"))}
      
      # get sensivity measures
      sensivity <- plyr::ddply(.data = errors, .variables = c("id"), .fun = get_sensitivity, rosetta_scores)
      rownames(sensivity) <- sensivity$id 
      sensivity$set <- set
      sensivity$grp <- grp

      
      sensivity <- sensivity[rnas, ]
      ifelse(is.null(sensivities), sensivities <- sensivity, sensivities <- rbind(sensivities, sensivity))
    }
  }
  return(sensivities)
}

combine_sensivities <- function(sensivities){
  # combines sensitivity data from scaha and rosetta
  sensivities$type <- "scaha"
  sensivities_scaha <- sensivities[, c("id", "roc", "nslr", "grp", "set")]
  sensivities_rosetta <- sensivities[, c("id","rosetta_roc", "nslr_rosetta", "grp", "set")]
  sensivities_rosetta$grp <- "rosetta"
  sensivities_rosetta <- unique(sensivities_rosetta)
  colnames(sensivities_rosetta) <- c("id", "roc", "nslr", "grp", "set")
  sensivities <- rbind(sensivities_scaha, sensivities_rosetta)
  
  sensivities$box_group <- 1
  sensivities$box_group[sensivities$set == "all" & sensivities$grp == "luc" ] <- 1
  sensivities$box_group[sensivities$set == "all" & sensivities$grp == "ruc" ] <- 2
  sensivities$box_group[sensivities$set == "all" & sensivities$grp == "rosetta" ] <- 3
  
  sensivities$box_group[sensivities$set == "training" & sensivities$grp == "luc" ] <- 4
  sensivities$box_group[sensivities$set == "training" & sensivities$grp == "ruc" ] <- 5
  sensivities$box_group[sensivities$set == "training" & sensivities$grp == "rosetta" ] <- 6
  
  sensivities$box_group[sensivities$set == "testing" & sensivities$grp == "luc" ] <- 7
  sensivities$box_group[sensivities$set == "testing" & sensivities$grp == "ruc" ] <- 8
  sensivities$box_group[sensivities$set == "testing" & sensivities$grp == "rosetta" ] <- 9
  sensivities$box_group <- as.factor(sensivities$box_group)
  
  return(sensivities)
}

merge_rosetta_scores <- function(rnas = unlist(strsplit("1ZC5 2LUN 2M4W 2MNC 2N6T 2M5U 2MXL 2NCI 1UUU 2JWV 2M12 2N2O 2QH2 1XHP 2K66 2LQZ 2M21 2MEQ 2N2P 2N7X 2QH4 2M22 2MFD 2N4L 1KKA 1PJY 2LI4 2LUB 4A4S 1L1W 2LBJ 2LK3 4A4T 1LC6 2FDT 2LBL 2LP9 2LV0 4A4U 1LDZ 2LDL 2LPA 1NC0 2LDT 1OW9 1YSV 2KOC 2LHP 2LU0 1SCL 1Z2J 2M24 2N6S", " "))){
  # merge all the rosetta scores for the RNAs in the challenge set into a single data.frame
  scores <- NULL
  rnas <- unique(rnas)
  for (rna in rnas){
    fname <- paste("rosetta_scores/", rna, ".txt", sep = "")
    try (silent = TRUE, {
      s <- read.table(fname, header = T)
      ifelse(is.null(scores), scores <- s, scores <- rbind(scores, s))
    })
    
  }
  return(scores)
}


make_sensitivity_boxplots <- function(s, fname = "Rplot.pdf", formula = nslr~box_group, scale = 1.0, width = 6.1, height = 2.6){
  # make sensitivity boxplots
  pdf(file = fname, width = scale*width, height = scale*height)
  {
    par(mfrow = c(1,3), cex = 0.8, mai = c(0.35, 0.35, 0.35, 0.35))
    col <- c("white", "grey", "grey25")
    p_all <- boxplot(formula, data = droplevels(subset(s, set=="all")), col = col, ylim = c(0, 1))
    abline(h=0.5, lwd = "2", lty = 2, col = "red")
    p_training <- boxplot(formula, data = droplevels(subset(s, set=="training")), col = col, ylim = c(0, 1))
    abline(h=0.5, lwd = "2", lty = 2, col = "red")
    p_testing <- boxplot(formula, data = droplevels(subset(s, set=="testing")), col = col, ylim = c(0, 1))
    abline(h=0.5, lwd = "2", lty = 2, col = "red")
  }
  dev.off()
  return(list(box_all=p_all, box_train = p_training, box_testing = p_testing))
}
