# 0 - goto working directory
setwd("~/Desktop/summer-papers/analysis/")

# 1 - source custom function
source("analysis_functions.R")

# 2 - load accuracy data
load("analysis_data.RData")

# 3 - get accuracy
acc <- combine_accuracy_data(data)

# 4 - accuracy
summarize_accuracy_for_paper(fname="proton_accuracy_raw.pdf", nuclei = "proton", codes = c("luu", "ruu"), lim = c(0, 1))
summarize_accuracy_for_paper(fname="carbon_accuracy_raw.pdf", nuclei = "carbon", codes = c("luu", "ruu"), lim = c(0, 4))
summarize_accuracy_for_paper(fname="carbon_accuracy_corrected.pdf", nuclei = "carbon", codes = c("luc", "ruc"), lim = c(0, 4))

# 5 - outlier inspection
# uncorrected proton outliers
# LARMORD: 2LUN (0.98 ppm), 2N2O (0.95 ppm), 2M22 (0.95 ppm), 2LV0 (0.49 ppm), 2M24 (0.30 ppm)
# RAMSEY: 2N7X (0.30 ppm), 2M24 (0.30 ppm)

# uncorrected carbon outliers
# LARMORD: 1Z2J (9.41 ppm), 2LV0 (8.46 ppm), 2LUN (6.53 ppm)
# RAMSEY: 2M24 (4.24 ppm)

# corrected carbon outliers
# LARMORD: 1Z2J (9.41 ppm 2LU0 (8.87 ppm), 2LV0 (8.46 ppm), 2LUN (6.53 ppm), 2N6T (6.26 ppm), and 2N2O(3.35 ppm),
# RAMSEY:  2N6T (4.30 ppm), 2LU0 (2.50 ppm), 1LC6 (2.11 ppm), and 1Z2J (1.68 ppm)

# 6a - do native-like structures have the lowest assignment errors  (actual vs. assigned vs)
actual_v_assigned_rmsd_1 <-get_rmsd_matrix(data, ref = "actual", comp = "assigned", N=1)
assigned_v_predicted_rmsd_1 <-get_rmsd_matrix(data, ref = "assigned", comp = "predicted", N=1)
actual_v_predicted_rmsd_1 <-get_rmsd_matrix(data, ref = "actual", comp = "predicted", N=1)

actual_v_assigned_rmsd_5 <-get_rmsd_matrix(data, ref = "actual", comp = "assigned", N=25)
assigned_v_predicted_rmsd_5 <-get_rmsd_matrix(data, ref = "assigned", comp = "predicted", N=25)
actual_v_predicted_rmsd_5 <-get_rmsd_matrix(data, ref = "actual", comp = "predicted", N=25)

# Levelplots
width <- 2.5
height <- 7.5
at <- seq(0, 7, 1)
pdf(file="actual_v_assigned.pdf", width = width, height = height)
make_levelplot(actual_v_assigned_rmsd_1, cols = c("white","black"), at = at)
dev.off()

pdf(file="assigned_v_predicted.pdf", width = width, height = height)
make_levelplot(assigned_v_predicted_rmsd_1, cols = c("white","black"), at = at)
dev.off()

pdf(file="actual_v_predicted.pdf", width = width, height = height)
make_levelplot(actual_v_predicted_rmsd_1, cols = c("white","black"), at = at)
dev.off()

# Fraction lower than threshold
summarize_identification_accuracy(matrices = c("actual_v_assigned_rmsd_1", "assigned_v_predicted_rmsd_1", "actual_v_predicted_rmsd_1"))
summarize_identification_accuracy(matrices = c("actual_v_assigned_rmsd_5", "assigned_v_predicted_rmsd_5", "actual_v_predicted_rmsd_5"))


# Get timings
timings <- make_timing_plot()
timings[,c("id", "bmrb", "conformers", "n_actual", "n_pred")]

# look at correlation between the RMSD and error
plot_rmsd_correlations(data, grp = "luc", corel_type = "rho", set = "all", col = "black", addplot = FALSE)
plot_rmsd_correlations(data, grp = "luc", corel_type = "rho", set = "training", col = "blue", addplot = TRUE)
plot_rmsd_correlations(data, grp = "luc", corel_type = "rho", set = "testing", col = "red", addplot = TRUE)

plot_rmsd_correlations(data, grp = "ruc", corel_type = "rho", set = "all", col = "black", addplot = FALSE)
plot_rmsd_correlations(data, grp = "ruc", corel_type = "rho", set = "training", col = "blue", addplot = TRUE)
plot_rmsd_correlations(data, grp = "ruc", corel_type = "rho", set = "testing", col = "red", addplot = TRUE)


# correlation distributions
save_correlation_boxplots(data, ref = "assigned", comp = "actual")

# assess sensitivity
s <- combine_sensivities(sensivities = get_sensivities(data, sets = c("all", "training", "testing"), grps = c("luc", "ruc"), rosetta_scores = TRUE))
make_sensitivity_boxplots(s, formula = nslr~box_group, scale = 1.0, fname = "sensitivity_nslr.pdf")
make_sensitivity_boxplots(s, formula = roc~box_group, scale = 1.0, fname = "sensitivity_auc.pdf")

# peek at data
pyshifts_ready_data(data, set = "luu", "1Z2J", cs="actual", filename = "1Z2J_observed.txt", predicted_file = FALSE)
pyshifts_ready_data(data, set = "luu", "1Z2J", cs="assigned", filename = "1Z2J_predicted.txt", predicted_file = TRUE)
