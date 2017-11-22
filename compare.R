lsap_assignment <- read.table("tests/assigned_shifts_1_1.txt")
custom_assignment <- read.table("tests_custom/assigned_shifts_1_1.txt")

cat("MAE between LSAP and Custom Assignments:\n")
cat(mean(abs(lsap_assignment$V5 - custom_assignment$V5)))
cat("\n")

cat("MAE between LSAP and Actual Assignments:\n")
cat(mean(abs(lsap_assignment$V5 - lsap_assignment$V6)))
cat("\n")

cat("MAE between Custom and Actual Assignments:\n")
cat(mean(abs(custom_assignment$V5 - custom_assignment$V6)))
cat("\n")