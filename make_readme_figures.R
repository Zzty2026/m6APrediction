# ===========================================================================
#  make_readme_figures.R
#
#  Creates the two images used in README.md:
#      man/figures/roc_curve.png
#      man/figures/prc_curve.png
#  and prints the AUC / PRC-AUC values to paste into README.md.
#
#  HOW TO RUN: open the m6APrediction project in RStudio, then run
#      source("make_readme_figures.R")
#  Run it from the PACKAGE ROOT (the folder that contains DESCRIPTION).
#  Do NOT put this script inside R/ - only the three functions live there.
#
#  Where do the numbers come from?
#  The random forest stores OUT-OF-BAG predictions for every training site
#  (rf$votes). OOB predictions are made by the trees that did NOT see that
#  site during training, so they are legitimate held-out predictions: no
#  leakage, and no need to re-load the Practical 4 test set.
#  If you already have the ROC/PRC plots from Practical 4, you may simply
#  copy them into man/figures/ instead of running this script.
# ===========================================================================

for (p in c("pROC", "PRROC")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

stopifnot(file.exists("DESCRIPTION"))   # we are in the package root
dir.create("man/figures", recursive = TRUE, showWarnings = FALSE)

ml_fit <- readRDS("inst/extdata/rf_fit.rds")

# Support both a plain randomForest object and a caret::train wrapper
rf <- if (!is.null(ml_fit$finalModel)) ml_fit$finalModel else ml_fit

if (is.null(rf$votes)) {
  stop("This model object has no $votes; use the Practical 4 plots instead, ",
       "or recompute predictions on a held-out set.")
}

truth <- rf$y                                  # true class labels
prob  <- as.numeric(rf$votes[, "Positive"])    # OOB probability of m6A

cat("Classes:", paste(levels(truth), collapse = " / "), "\n")
cat("Sites:", length(truth),
    "| positives:", sum(truth == "Positive"), "\n")

# ---- ROC ------------------------------------------------------------------
roc_obj <- pROC::roc(response = truth, predictor = prob,
                     levels = c("Negative", "Positive"), quiet = TRUE)
auc_val <- as.numeric(pROC::auc(roc_obj))

png("man/figures/roc_curve.png", width = 1600, height = 1300, res = 220)
par(mar = c(4.5, 4.5, 2, 1))
plot(roc_obj, col = "#1f6fb2", lwd = 2.5, legacy.axes = TRUE,
     print.auc = TRUE, print.auc.cex = 1.1, print.auc.y = 0.35,
     auc.polygon = TRUE, auc.polygon.col = "#dbe8f5",
     grid = c(0.1, 0.2), grid.col = "grey85",
     xlab = "False positive rate (1 - specificity)",
     ylab = "True positive rate (sensitivity)",
     main = "ROC curve - m6A random forest (out-of-bag)")
dev.off()

# ---- Precision-recall -----------------------------------------------------
prc_obj <- PRROC::pr.curve(scores.class0 = prob[truth == "Positive"],
                           scores.class1 = prob[truth == "Negative"],
                           curve = TRUE)
prc_auc <- prc_obj$auc.integral

png("man/figures/prc_curve.png", width = 1600, height = 1300, res = 220)
par(mar = c(4.5, 4.5, 2, 1))
plot(prc_obj, color = "#1f6fb2", lwd = 2.5,
     xlab = "Recall (sensitivity)",
     ylab = "Precision",
     main = "Precision-recall curve - m6A random forest (out-of-bag)")
abline(h = mean(truth == "Positive"), lty = 2, col = "grey40")
legend("bottomleft",
       legend = c("PR curve", "Prevalence baseline"),
       col = c("#1f6fb2", "grey40"), lty = c(1, 2), lwd = c(2.5, 1),
       bty = "n", cex = 0.9)
dev.off()

# ---- Numbers for README.md ------------------------------------------------
cat("\n--- paste these into README.md -------------------------------\n")
cat(sprintf("AUC = %.2f ; PRC-AUC = %.2f\n", auc_val, prc_auc))
cat("-------------------------------------------------------------\n")
cat("Files written: man/figures/roc_curve.png, man/figures/prc_curve.png\n")
