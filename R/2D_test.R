library(nmR)

setwd("~/GitSoftware/SCAHA/")
names <- c("model", "resid", "resname", "nucleus", "predCS", "expCS")
file <- "tests/larmord_2KOC_single.txt"
tmp <- nmR::load_cs_data(csfile=file, accuracyFile = "tests/larmord_accuracy_resname_nucleus.txt", atomBasedWeights = TRUE, names = names)

larmord_peaks <- nmR::create_peaks(tmp)
hmqc_peaks <- read.table(file = "tests/region5_2D_hmqc.tab", header = TRUE)
hmqc_peaks <- hmqc_peaks[hmqc_peaks$volume>0 & hmqc_peaks$volume<30, ]



