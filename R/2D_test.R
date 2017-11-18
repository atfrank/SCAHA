setwd("~/GitSoftware/SCAHA/")

library(nmR)
source("R/library.R")

names <- c("model", "resid", "resname", "nucleus", "predCS", "expCS")
file <- "tests/region5_2D_larmord.dat"
tmp <- nmR::load_cs_data(csfile=file, accuracyFile = "tests/larmord_accuracy_resname_nucleus.txt", atomBasedWeights = TRUE, names = names)

# loads assigned computed chemical shift data
names <- c("conformation", "resid", "type_H", "type_C", "peak_H", "peak_C", "weight_H", "weight_C")
assgined_computed_cs <- nmR::create_peaks(tmp)
colnames(assgined_computed_cs) <- names

# load unassigned chemical shift data
unassgined_data <- read.table(file = "tests/region5_2D_hmqc.tab", header = TRUE)
unassgined_data <- unassgined_data[unassgined_data$volume>0, ]
unassgined_data <- duplicate_unassigned_data(unassgined_data, assgined_computed_cs)

costmat <- matrix(0, nrow = nrow(unassgined_data), ncol = nrow(assgined_computed_cs))
for (i in 1:nrow(unassgined_data)){
  d1 <- matrix(as.vector(unlist(unassgined_data[i, c("ppm1", "ppm2")])), ncol=2, nrow=nrow(assgined_computed_cs), byrow = TRUE)
  d2 <- matrix(as.vector(unlist(assgined_computed_cs[, c("peak_H", "peak_C")])), ncol=2, byrow = FALSE)
  w2 <- matrix(as.vector(unlist(assgined_computed_cs[, c("weight_H", "weight_C")])), ncol=2, byrow = FALSE)
  w1 <- unassgined_data$volume[i]*unassgined_data$intensity[i]
  #costmat[i, ] <- rowSums(w2*w2*(d1-d2)*(d1-d2))
  costmat[i, ] <- (1/w1)*(w2[,1]*abs(d1[,1]-d2[,1])+w2[,2]*abs(d1[,2]-d2[,2]))
  #costmat[i, ] <- abs(d1[,1]-d2[,1])+abs(d1[,2]-d2[,2])
}

a <- solve_LSAP(t(costmat))
assgined_computed_cs$assigned_peak_H <-unassgined_data$ppm1[a]
assgined_computed_cs$assigned_peak_C <-unassgined_data$ppm2[a]

