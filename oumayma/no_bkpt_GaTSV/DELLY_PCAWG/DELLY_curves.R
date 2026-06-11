# -------------------------------
# Load libraries
# -------------------------------
library(e1071)
library(caret)
library(pROC)
library(ROCR)
library(data.table)
library(tibble)
library(yardstick)
library(dplyr)
set.seed(42)

# -------------------------------
# Load data & Prepare labels
# -------------------------------
# Load unscaled combined data to get truth labels
unscaled_data <- readRDS("/data/beroukhim1/oumayma/no_bkpt_GaTSV/DELLY_PCAWG/output/processed/combined_delly.rds")
unscaled_data <- as.data.table(unscaled_data)
unscaled_data <- unscaled_data[SPAN >= 1e3 | SPAN == -1, ]

# Load scaled features
df <- readRDS("/data/beroukhim1/oumayma/no_bkpt_GaTSV/DELLY_PCAWG/output/delly_test_scaled.rds")
df <- as.data.table(df)

# Attach truth labels
df$sv_class <- ifelse(unscaled_data$truth == "SOMATIC", 1, 0)

# -------------------------------
# Define test set indices 
# -------------------------------
y_full <- factor(df$sv_class, levels = c(1,0), labels = c("somatic", "germline"))
trainIndex <- createDataPartition(y_full, p = 0.66, list = FALSE)
y_test <- y_full[-trainIndex]

# -------------------------------
# Prepare test matrix 
# -------------------------------
exclude_features <- c("log_homlen", "log_insertion_len", "hom_gc", "insertion_gc", "sv_class")
X_test <- as.data.frame(df[-trainIndex, setdiff(names(df), exclude_features), with = FALSE])

# -------------------------------
# Load pretrained SVM model safely
# -------------------------------
tmp_env <- new.env()
load("/data/beroukhim1/oumayma/no_bkpt_GaTSV/GaTS_BpA_031826.rda", envir = tmp_env)
GaTSV <- tmp_env[[ls(tmp_env)[1]]]
rm(tmp_env)

# -------------------------------
# Predict probabilities and classes
# -------------------------------
cutoff <- 0.2684

# Get probabilities
probs_matrix <- attr(predict(GaTSV, X_test, probability = TRUE), "probabilities")

# Safely extract the positive class ("1" or "SOMATIC")
pos_col <- grep("somatic|1", colnames(probs_matrix), value = TRUE, ignore.case = TRUE)
if(length(pos_col) == 0) pos_col <- colnames(probs_matrix)[2] else pos_col <- pos_col[1]

probs <- as.numeric(probs_matrix[, pos_col])

# Get class predictions based on cutoff
preds <- factor(ifelse(probs > cutoff, "somatic", "germline"), levels = c("somatic","germline"))

# -------------------------------
# Calculate metrics
# -------------------------------
calc_stats <- function(truth, preds, probs) {
  tpr  <- sens_vec(truth, preds, event_level = "first")
  spec <- spec_vec(truth, preds, event_level = "first")
  ppv  <- ppv_vec(truth, preds, event_level = "first")
  auc  <- roc_auc_vec(truth, probs, event_level = "first")
  tibble(AUC = auc, PPV = ppv, TPR = tpr, Specificity = spec)
}

results <- calc_stats(y_test, preds, probs) %>% mutate(Model = "GaTSV (Delly)")

cat("\n--- Test Set Results ---\n")
print(results)


# ======================================================
#        ROC + PR Curves (Saved to PDF)
# ======================================================
pdf("/data/beroukhim1/oumayma/no_bkpt_GaTSV/DELLY_PCAWG/output/delly_svm_curves.pdf", width = 10, height = 5)

par(mfrow = c(1, 2),
    mar = c(5, 5, 4, 2),
    oma = c(0, 0, 2, 0),
    pty = "s")

# -------------------------------
# 1. ROC Curve
# -------------------------------
roc_obj <- roc(y_test, probs, levels = c("germline", "somatic"))

plot(roc_obj, col = "steelblue", lwd = 3, xlim = c(1, 0),
     main = "SVM on Delly Test-Set: ROC Curve", asp = 1)

# AUC text (bottom-left)
text(0.8, 0.15, sprintf("GaTSV AUC = %.3f", results$AUC),
     col = "steelblue", cex = 1.2, font = 2, adj = 0)

# -------------------------------
# 2. PR Curve
# -------------------------------
pred_obj <- ROCR::prediction(probs, as.numeric(y_test == "somatic"))
perf_pr  <- ROCR::performance(pred_obj, "prec", "rec")
pr_auc   <- ROCR::performance(pred_obj, "aucpr")@y.values[[1]]

# Plot PR curve
plot(perf_pr@x.values[[1]], perf_pr@y.values[[1]],
     type = "l", col = "darkorange", lwd = 3,
     main = "SVM on Delly Test-Set: PR Curve",
     xlab = "Recall", ylab = "Precision",
     xlim = c(0, 1), ylim = c(0, 1), asp = 1)

# Add baseline (random guess line)
abline(h = mean(as.numeric(y_test == "somatic")), lty = 2, col = "gray")

# PR AUC text (bottom-left)
text(0.1, 0.15, sprintf("GaTSV PR AUC = %.3f", pr_auc),
     col = "darkorange", cex = 1.2, font = 2, adj = 0)

# Reset layout and close PDF device
par(mfrow = c(1, 1))
dev.off()

cat("\nCurves successfully saved to: /data/beroukhim1/oumayma/no_bkpt_GaTSV/DELLY_PCAWG/output/delly_svm_curves.pdf\n")