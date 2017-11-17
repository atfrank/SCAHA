library(nmR)
create_peaks <- function(cs, protons=c("H1'", "H2'", "H3'", "H4'", "H5'", "H5''", "H2", "H5", "H6", "H8"), carbons=c("C1'", "C2'", "C3'", "C4'", "C5'", "C5'", "C2", "C5", "C6", "C8")){
  peak_H <- NULL
  peak_C <- NULL
  type_H <- NULL
  type_C <- NULL
  weight_H <- NULL
  weight_C <- NULL

  for (i in 1:length(protons)){
    if(length(cs$predCS[cs$nucleus==protons[i]])!=0){
      peak_H <- c(peak_H, cs$predCS[cs$nucleus==protons[i]])
      peak_C <- c(peak_C, cs$predCS[cs$nucleus==carbons[i]])
      type_H <- c(type_H, protons[i])
      type_C <- c(type_C, carbons[i])
      weight_H <- c(weight_H, cs$weight[cs$nucleus==protons[i]])
      weight_C <- c(weight_C, cs$weight[cs$nucleus==carbons[i]])
    }
  }
  return(data.frame(type_H=type_H, type_C=type_C, peak_H=peak_H, peak_C=peak_C, weight_H=weight_H, weight_C=weight_C))
}


cstable2peaks <- function(cs, grouping=c("model", "resid")){
  plyr::ddply(.data = cs, .variables = grouping, .fun = create_peaks)
}


names <- c("model", "resid", "resname", "nucleus", "predCS", "expCS")
file <- "tests/larmord_2KOC_single.txt"
tmp <- nmR::load_cs_data(csfile=file, accuracyFile = "tests/larmord_accuracy_resname_nucleus.txt", atomBasedWeights = TRUE, names = names)

cstable2peaks(tmp)