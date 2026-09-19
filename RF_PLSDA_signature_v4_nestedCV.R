# Title: [RF_PLSDA_signature.R]
# Author: [Murli Nair]
# Date: [2026-03-15]
# License: GNU General Public License v3.0+
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.
#
# USER RESPONSIBILITY:
# You are solely responsible for validating, testing, and verifying
# this code for your specific use case. See LICENSE.md for full terms.

# For reproducibility use the version of the packages provide in SessionInfo.txt

############################################################
## RF + PLS-DA + Venn + Heatmap of discriminant signature ##
############################################################

## 0. Packages ----
pkgs <- c(
  "dplyr", "tibble", "matrixStats",
  "randomForest", "mixOmics",
  "VennDiagram", "pheatmap", "pROC"
)
# missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
# if (length(missing)) {
#   install.packages(setdiff(missing, c("mixOmics","randomForest","VennDiagram")))
#   if ("mixOmics" %in% missing || "randomForest" %in% missing || "VennDiagram" %in% missing) {
#     if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#     BiocManager::install(setdiff(missing, pkgs[!pkgs %in% c("mixOmics","randomForest","VennDiagram")]), ask = FALSE)
#   }
# }
invisible(lapply(pkgs, library, character.only = TRUE))

set.seed(1234)

############################################################
## 1. LOAD / PREPARE DATA                                ##
############################################################
# First column = SYMBOL, rest = samples
# expr_df <- mergedData %>%
#   dplyr::select(SYMBOL, MDAMB157_CON1:APCshRNA2_PTX3)

expr_df <- read.csv("normalizedCounts_withSymbols_andENSEMBL2.csv")
expr_df$SYMBOL[duplicated(expr_df$SYMBOL)]<-paste(expr_df$SYMBOL[duplicated(expr_df$SYMBOL)],
                                                  expr_df$ENSEMBL[which(duplicated(expr_df$SYMBOL)==TRUE)], sep="-")

expr_df<-expr_df[,c(-2,(length(expr_df[1,])*-1))]
head(expr_df)
# Extract sample names (all except gene column)
gene_col_name <- colnames(expr_df)[1]
samples <- colnames(expr_df)[-1]


# Build annotation: Genotype + Treatment + Group (Genotype_Treatment)
Genotype  <- gsub("_.*", "", samples)
Treatment <- gsub(".*_(CON|CIS|PTX)\\d+$", "\\1", samples)

ann <- data.frame(
  Genotype  = Genotype,
  Treatment = Treatment,
  Group     = paste(Genotype, Treatment, sep = "_"),
  row.names = samples
)

############################################################
## 2. BUILD EXPRESSION MATRIX (genes x samples)          ##
############################################################

# Coerce to matrix: rows = genes, cols = samples
expr_mat <- expr_df[, -1, drop = FALSE]
rownames(expr_mat) <- expr_df[[1]]
expr_mat <- as.matrix(expr_mat)
mode(expr_mat) <- "numeric"

# Optional: log1p-transform (if counts/TPM-like)
#expr_mat_log <- log1p(expr_mat)
expr_mat_log <- log2(expr_mat + 1)
# Filter low-variance genes to stabilize models
vars <- matrixStats::rowVars(expr_mat_log, na.rm = TRUE)
var_thresh <- quantile(vars, 0.25)  # keep top 75% variable genes; adjust as needed
keep <- vars > var_thresh
expr_filt <- expr_mat_log[keep, , drop = FALSE]

cat("Genes before filtering:", nrow(expr_mat_log), "\n")
cat("Genes after filtering :", nrow(expr_filt), "\n")

############################################################
## 3. RANDOM FOREST FEATURE IMPORTANCE                  ##
############################################################

# Random Forest expects: samples as rows, genes as columns
X_rf <- t(expr_filt)

# Response: here we classify Treatment (CON/CIS/PTX)
y_rf <- factor(ann[colnames(expr_filt), "Treatment"])

rf_model <- randomForest::randomForest(
  x = X_rf,
  y = y_rf,
  importance = TRUE,
  ntree = 2000
)

# Importance: Mean Decrease Accuracy (type = 1)
rf_imp <- randomForest::importance(rf_model, type = 1)
rf_imp <- as.data.frame(rf_imp)
rf_imp$gene <- rownames(rf_imp)

# Top N genes by RF importance
topN <- 100  # you can tune this
rf_top <- rf_imp %>%
  arrange(desc(MeanDecreaseAccuracy)) %>%
  slice(1:topN)

rf_genes <- rf_top$gene

cat("Top RF genes:", length(rf_genes), "\n")

############################################################
## 4. PLS-DA + VIP SCORES (mixOmics)                    ##
############################################################

# mixOmics expects: samples x variables
X_pls <- X_rf  # same as RF input
Y_pls <- y_rf

# Choose number of components (2–3 usually enough)
ncomp <- 2
pls_model <- mixOmics::plsda(X_pls, Y_pls, ncomp = ncomp)

# VIP scores for all genes across components
vip_mat <- mixOmics::vip(pls_model)  # nvar x ncomp
vip_df <- as.data.frame(vip_mat)
vip_df$gene <- rownames(vip_df)

# Summarize VIP across components (e.g., mean VIP)
vip_df$VIP_mean <- rowMeans(vip_df[, 1:ncomp, drop = FALSE])

# Top N genes by VIP
pls_top <- vip_df %>%
  arrange(desc(VIP_mean)) %>%
  slice(1:topN)

pls_genes <- pls_top$gene

cat("Top PLS-DA genes:", length(pls_genes), "\n")

############################################################
## 5. OVERLAP (SIGNATURE GENES) + VENN DIAGRAM          ##
############################################################

signature_genes <- intersect(rf_genes, pls_genes)
cat("Signature genes (RF ∩ PLS-DA):", length(signature_genes), "\n")

# Venn diagram: RF vs PLS
venn_list <- list(
  RF     = rf_genes,
  PLS_DA = pls_genes
)

# Draw to current device OR save to file
grid.newpage()
venn_plot <- VennDiagram::venn.diagram(
  x = venn_list,
  filename = NULL,
  fill = c("#377EB8", "#E41A1C"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.3,
  main = "Overlap of RF and PLS-DA top genes"
)
grid::grid.draw(venn_plot)

# Optionally save venn to file:
# Save it using grDevices::dev.copy()
grDevices::dev.copy(pdf, file = "RF_PLSDA_venn_v3.pdf", width = 10, height = 10)
grDevices::dev.off()
while (!is.null(dev.list())) dev.off()

############################################################
## 6. HEATMAP OF SIGNATURE GENES                        ##
############################################################

if (length(signature_genes) > 1) {
  expr_sig <- expr_mat_log[signature_genes, , drop = FALSE]
  
  # z-score per gene for heatmap
  expr_sig_z <- t(scale(t(expr_sig)))
  expr_sig_z[!is.finite(expr_sig_z)] <- 0
  
  # annotation for columns (samples)
  ann_col <- ann[colnames(expr_sig_z), c("Genotype","Treatment"), drop = FALSE]
  
  # Heatmap
  pheatmap::pheatmap(
    expr_sig_z,
    annotation_col = ann_col,
    show_rownames = TRUE,
    show_colnames = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    main = sprintf("Signature genes (RF ∩ PLS-DA), n = %d", length(signature_genes)),
    fontsize_row = 6,
    fontsize_col = 8
  )
  
  # Optionally save to file:
  pheatmap::pheatmap(expr_sig_z,
                     annotation_col = ann_col,
                     show_rownames = TRUE,
                     show_colnames = TRUE,
                     clustering_distance_rows = "euclidean",
                     clustering_distance_cols = "euclidean",
                     clustering_method = "complete",
                     main = sprintf("Signature genes (RF & PLS-DA), n = %d", length(signature_genes)),
                     fontsize_row = 6,
                     fontsize_col = 8,
                     filename = "signature_heatmap_v3.pdf",  # ← Uses correct device + handles dev.off()
                     width = 8,   # inches
                     height = 10  # inches
                     )
  
  
  
  
} else {
  cat("Signature gene set is too small for a meaningful heatmap.\n")
}


############################################################
## 6b. HEATMAP OF SIGNATURE GENES - NORMALIZED EXPRESSION  ##
############################################################

if (length(signature_genes) > 1) {
  
  # log1p-transformed normalized counts
  expr_sig <- expr_mat_log[
    signature_genes,
    ,
    drop = FALSE
  ]
  
  # Column annotation
  ann_col <- ann[
    colnames(expr_sig),
    c("Genotype", "Treatment"),
    drop = FALSE
  ]
  
  # --------------------------------------------------------
  # Order samples by GENOTYPE first, then TREATMENT
  # --------------------------------------------------------
  
  genotype_order <- c(
    "MDAMB157",
    "APCshRNA1",
    "APCshRNA2"
  )
  
  treatment_order <- c(
    "CON",
    "CIS",
    "PTX"
  )
  
  sample_order <- order(
    factor(
      ann_col$Genotype,
      levels = genotype_order
    ),
    factor(
      ann_col$Treatment,
      levels = treatment_order
    )
  )
  
  expr_sig <- expr_sig[
    ,
    sample_order,
    drop = FALSE
  ]
  
  ann_col <- ann_col[
    sample_order,
    ,
    drop = FALSE
  ]
  
  # --------------------------------------------------------
  # Heatmap of actual log-transformed normalized counts
  #
  # NO gene-wise z-score
  # NO column clustering
  # --------------------------------------------------------
  
  pheatmap::pheatmap(
    expr_sig,
    
    annotation_col = ann_col,
    
    show_rownames = TRUE,
    show_colnames = TRUE,
    
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    
    clustering_distance_rows = "euclidean",
    clustering_method = "complete",
    
    main = sprintf(
      "RF/PLS-DA treatment-discriminant genes, n = %d",
      length(signature_genes)
    ),
    
    fontsize_row = 6,
    fontsize_col = 7,
    
    angle_col = 45,
    
    filename =
      "signature_heatmap_normalized_counts_v3.pdf",
    
    width = 10,
    height = 10
  )
}






############################################################
## 7. EXPORT SIGNATURE TABLE                            ##
############################################################

signature_table <- expr_df %>%
  dplyr::filter(.data[[gene_col_name]] %in% signature_genes)

# Add RF + VIP scores for context
signature_table <- signature_table %>%
  left_join(rf_top[, c("gene","MeanDecreaseAccuracy")],
            by = c(SYMBOL = "gene")) %>%
  left_join(pls_top[, c("gene","VIP_mean")],
            by = c(SYMBOL = "gene"))
#Save the signature table as a file
 write.csv(signature_table, "signature_genes_RF_PLSDA_v3.csv", row.names = FALSE)

cat("Pipeline complete.\n")

############################################################
## 8. LEAKAGE-SAFE REPEATED STRATIFIED CROSS-VALIDATION  ##
############################################################
# Purpose:
#   Evaluate whether the RF/PLS-DA treatment-discriminant pipeline
#   generalizes to held-out samples while preventing information leakage.
#
# IMPORTANT:
#   Every preprocessing and feature-selection step below is performed
#   INSIDE each training fold. The held-out test samples are not used for
#   variance filtering, RF ranking, PLS-DA ranking, or model fitting.
#
# This validation is separate from the full-dataset 43-gene feature set
# generated above, which is retained for biological interpretation.

set.seed(20260916)

CV_K <- 3
CV_REPEATS <- 20
CV_TOPN <- 100
CV_NCOMP <- 2
CV_NTREE_FS <- 2000
CV_NTREE_FINAL <- 2000

# Permutation testing can be computationally expensive.
# N_PERM = 100 gives a minimum attainable empirical P value of ~0.0099.
# Increase to 1000 for a final analysis if computationally feasible.
N_PERM <- 1000
PERM_CV_REPEATS <- 1

# ----------------------------------------------------------
# 8a. Helper functions
# ----------------------------------------------------------

make_stratified_folds <- function(y, k = 3, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  y <- factor(y)
  folds <- vector("list", k)
  for (lev in levels(y)) {
    idx <- which(y == lev)
    idx <- sample(idx, length(idx), replace = FALSE)
    split_idx <- split(idx, rep(seq_len(k), length.out = length(idx)))
    for (j in seq_len(k)) {
      folds[[j]] <- c(folds[[j]], split_idx[[j]])
    }
  }
  lapply(folds, sort)
}

select_rf_pls_features <- function(X_train, y_train,
                                   topN = 100,
                                   ncomp = 2,
                                   ntree = 2000) {
  # X_train: samples x genes, already log2(normalized count + 1)

  # Variance filtering is calculated using TRAINING samples only.
  train_vars <- matrixStats::colVars(X_train, na.rm = TRUE)
  var_thresh <- quantile(train_vars, 0.25, na.rm = TRUE)
  keep_genes <- names(train_vars)[is.finite(train_vars) & train_vars > var_thresh]

  Xf <- X_train[, keep_genes, drop = FALSE]

  # Random Forest ranking on training data only.
  rf_fs <- randomForest::randomForest(
    x = Xf,
    y = factor(y_train),
    importance = TRUE,
    ntree = ntree
  )
  rf_imp <- randomForest::importance(rf_fs, type = 1)
  rf_imp <- data.frame(
    gene = rownames(rf_imp),
    MeanDecreaseAccuracy = rf_imp[, "MeanDecreaseAccuracy"],
    row.names = NULL
  )
  rf_imp <- rf_imp[order(rf_imp$MeanDecreaseAccuracy, decreasing = TRUE), ]
  rf_top_genes <- head(rf_imp$gene, min(topN, nrow(rf_imp)))

  # PLS-DA ranking on training data only.
  pls_fs <- mixOmics::plsda(Xf, factor(y_train), ncomp = ncomp)
  vip_mat <- mixOmics::vip(pls_fs)
  vip_mean <- rowMeans(vip_mat[, seq_len(ncomp), drop = FALSE])
  vip_df <- data.frame(
    gene = names(vip_mean),
    VIP_mean = as.numeric(vip_mean),
    row.names = NULL
  )
  vip_df <- vip_df[order(vip_df$VIP_mean, decreasing = TRUE), ]
  pls_top_genes <- head(vip_df$gene, min(topN, nrow(vip_df)))

  selected <- intersect(rf_top_genes, pls_top_genes)

  list(
    genes = selected,
    rf_top = rf_top_genes,
    pls_top = pls_top_genes,
    n_after_variance_filter = ncol(Xf)
  )
}

calc_multiclass_metrics <- function(obs, pred, prob, class_levels) {
  obs <- factor(obs, levels = class_levels)
  pred <- factor(pred, levels = class_levels)

  cm <- table(obs, pred)
  accuracy <- sum(diag(cm)) / sum(cm)

  sensitivity <- specificity <- setNames(rep(NA_real_, length(class_levels)), class_levels)

  for (cl in class_levels) {
    tp <- cm[cl, cl]
    fn <- sum(cm[cl, ]) - tp
    fp <- sum(cm[, cl]) - tp
    tn <- sum(cm) - tp - fn - fp

    sensitivity[cl] <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
    specificity[cl] <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
  }

  balanced_accuracy <- mean(sensitivity, na.rm = TRUE)

  # One-vs-rest AUC for each treatment class; macro-average across classes.
  auc_by_class <- setNames(rep(NA_real_, length(class_levels)), class_levels)
  for (cl in class_levels) {
    truth_bin <- as.integer(obs == cl)
    if (length(unique(truth_bin)) == 2) {
      roc_obj <- suppressMessages(
        pROC::roc(response = truth_bin,
                  predictor = prob[, cl],
                  levels = c(0, 1),
                  direction = "<",
                  quiet = TRUE)
      )
      auc_by_class[cl] <- as.numeric(pROC::auc(roc_obj))
    }
  }
  macro_auc <- mean(auc_by_class, na.rm = TRUE)

  list(
    accuracy = accuracy,
    balanced_accuracy = balanced_accuracy,
    sensitivity = sensitivity,
    specificity = specificity,
    auc_by_class = auc_by_class,
    macro_auc = macro_auc,
    confusion_matrix = cm
  )
}

run_repeated_cv <- function(X, y,
                            k = 3,
                            repeats = 20,
                            topN = 100,
                            ncomp = 2,
                            ntree_fs = 2000,
                            ntree_final = 2000,
                            seed = 1,
                            collect_features = TRUE,
                            verbose = TRUE) {

  y <- factor(y)
  class_levels <- levels(y)
  sample_ids <- rownames(X)

  repeat_metrics <- vector("list", repeats)
  all_predictions <- vector("list", repeats)
  selected_features_all <- list()
  fold_counter <- 0

  for (r in seq_len(repeats)) {
    if (verbose) cat("CV repeat", r, "of", repeats, "\n")

    folds <- make_stratified_folds(y, k = k, seed = seed + r)

    obs_r <- pred_r <- character(0)
    sample_r <- character(0)
    prob_r <- matrix(numeric(0), nrow = 0, ncol = length(class_levels),
                     dimnames = list(NULL, class_levels))

    for (f in seq_len(k)) {
      fold_counter <- fold_counter + 1
      test_idx <- folds[[f]]
      train_idx <- setdiff(seq_len(nrow(X)), test_idx)

      X_train <- X[train_idx, , drop = FALSE]
      X_test <- X[test_idx, , drop = FALSE]
      y_train <- droplevels(y[train_idx])
      y_test <- y[test_idx]

      fs <- select_rf_pls_features(
        X_train = X_train,
        y_train = y_train,
        topN = topN,
        ncomp = ncomp,
        ntree = ntree_fs
      )

      sel <- fs$genes
      if (length(sel) < 2) {
        warning(sprintf("Repeat %d fold %d selected fewer than 2 genes; fold skipped.", r, f))
        next
      }

      if (collect_features) {
        selected_features_all[[paste0("R", r, "_F", f)]] <- sel
      }

      # Final classifier trained ONLY on selected training-fold genes.
      rf_final <- randomForest::randomForest(
        x = X_train[, sel, drop = FALSE],
        y = y_train,
        ntree = ntree_final,
        importance = TRUE
      )

      pred <- predict(rf_final, X_test[, sel, drop = FALSE], type = "response")
      prob <- predict(rf_final, X_test[, sel, drop = FALSE], type = "prob")
      prob <- prob[, class_levels, drop = FALSE]

      obs_r <- c(obs_r, as.character(y_test))
      pred_r <- c(pred_r, as.character(pred))
      sample_r <- c(sample_r, sample_ids[test_idx])
      prob_r <- rbind(prob_r, prob)
    }

    metrics <- calc_multiclass_metrics(
      obs = obs_r,
      pred = pred_r,
      prob = prob_r,
      class_levels = class_levels
    )

    repeat_metrics[[r]] <- data.frame(
      cv_repeat = r,
      accuracy = metrics$accuracy,
      balanced_accuracy = metrics$balanced_accuracy,
      macro_auc = metrics$macro_auc,
      sensitivity_CON = metrics$sensitivity["CON"],
      sensitivity_CIS = metrics$sensitivity["CIS"],
      sensitivity_PTX = metrics$sensitivity["PTX"],
      specificity_CON = metrics$specificity["CON"],
      specificity_CIS = metrics$specificity["CIS"],
      specificity_PTX = metrics$specificity["PTX"],
      auc_CON = metrics$auc_by_class["CON"],
      auc_CIS = metrics$auc_by_class["CIS"],
      auc_PTX = metrics$auc_by_class["PTX"],
      stringsAsFactors = FALSE
    )
    
    pred_df <- data.frame(
      cv_repeat = r,
      sample = sample_r,
      observed = obs_r,
      predicted = pred_r,
      stringsAsFactors = FALSE
    )
    pred_df <- cbind(pred_df, as.data.frame(prob_r))
    all_predictions[[r]] <- pred_df
  }

  metrics_df <- do.call(rbind, repeat_metrics)
  predictions_df <- do.call(rbind, all_predictions)

  feature_frequency <- NULL
  if (collect_features && length(selected_features_all) > 0) {
    all_sel <- unlist(selected_features_all, use.names = FALSE)
    freq <- sort(table(all_sel), decreasing = TRUE)
    feature_frequency <- data.frame(
      gene = names(freq),
      selected_folds = as.integer(freq),
      selection_frequency = as.integer(freq) / length(selected_features_all),
      row.names = NULL
    )
  }

  list(
    metrics = metrics_df,
    predictions = predictions_df,
    feature_frequency = feature_frequency,
    selected_features_by_fold = selected_features_all
  )
}

summarize_metric_ci <- function(x) {
  x <- x[is.finite(x)]
  c(
    mean = mean(x),
    median = median(x),
    lower_95 = unname(quantile(x, 0.025, type = 8)),
    upper_95 = unname(quantile(x, 0.975, type = 8))
  )
}

# ----------------------------------------------------------
# 8b. Build sample x gene matrix for validation
# ----------------------------------------------------------
# Use the same normalized-count source as the exploratory analysis.
# The log transformation itself is unsupervised; variance filtering and
# all supervised feature selection are re-estimated inside each train fold.

X_cv <- t(log2(expr_mat + 1))
y_cv <- factor(ann[rownames(X_cv), "Treatment"], levels = c("CON", "CIS", "PTX"))

cat("\nValidation sample count:", nrow(X_cv), "\n")
cat("Class distribution:\n")
print(table(y_cv))

cv_results <- run_repeated_cv(
  X = X_cv,
  y = y_cv,
  k = CV_K,
  repeats = CV_REPEATS,
  topN = CV_TOPN,
  ncomp = CV_NCOMP,
  ntree_fs = CV_NTREE_FS,
  ntree_final = CV_NTREE_FINAL,
  seed = 20260916,
  collect_features = TRUE,
  verbose = TRUE
)

write.csv(cv_results$metrics,
          "RF_PLSDA_repeatedCV_metrics_v4.csv",
          row.names = FALSE)
write.csv(cv_results$predictions,
          "RF_PLSDA_repeatedCV_predictions_v4.csv",
          row.names = FALSE)
if (!is.null(cv_results$feature_frequency)) {
  write.csv(cv_results$feature_frequency,
            "RF_PLSDA_feature_selection_stability_v4.csv",
            row.names = FALSE)
}

# Summary with empirical 95% intervals across CV repeats.
metric_cols <- c(
  "accuracy", "balanced_accuracy", "macro_auc",
  "sensitivity_CON", "sensitivity_CIS", "sensitivity_PTX",
  "specificity_CON", "specificity_CIS", "specificity_PTX",
  "auc_CON", "auc_CIS", "auc_PTX"
)

cv_summary <- do.call(rbind, lapply(metric_cols, function(m) {
  vals <- summarize_metric_ci(cv_results$metrics[[m]])
  data.frame(
    metric = m,
    mean = vals["mean"],
    median = vals["median"],
    lower_95 = vals["lower_95"],
    upper_95 = vals["upper_95"],
    row.names = NULL
  )
}))

write.csv(cv_summary,
          "RF_PLSDA_repeatedCV_summary_v4.csv",
          row.names = FALSE)

cat("\nRepeated CV summary:\n")
print(cv_summary)

# ----------------------------------------------------------
# 8c. Permutation test
# ----------------------------------------------------------
# Treatment labels are permuted and the COMPLETE CV pipeline is rerun,
# including within-fold variance filtering and RF/PLS-DA feature selection.
# We use mean balanced accuracy across repeats as the permutation statistic.

observed_stat <- mean(cv_results$metrics$balanced_accuracy, na.rm = TRUE)
perm_stats <- rep(NA_real_, N_PERM)

set.seed(20260917)
for (b in seq_len(N_PERM)) {
  cat("Permutation", b, "of", N_PERM, "\n")
  y_perm <- sample(y_cv, length(y_cv), replace = FALSE)
  y_perm <- factor(y_perm, levels = levels(y_cv))

  perm_res <- run_repeated_cv(
    X = X_cv,
    y = y_perm,
    k = CV_K,
    repeats = PERM_CV_REPEATS,
    topN = CV_TOPN,
    ncomp = CV_NCOMP,
    ntree_fs = CV_NTREE_FS,
    ntree_final = CV_NTREE_FINAL,
    seed = 910000 + b,
    collect_features = FALSE,
    verbose = FALSE
  )

  perm_stats[b] <- mean(perm_res$metrics$balanced_accuracy, na.rm = TRUE)
}

perm_p <- (1 + sum(perm_stats >= observed_stat, na.rm = TRUE)) /
          (1 + sum(is.finite(perm_stats)))

perm_out <- data.frame(
  permutation = seq_len(N_PERM),
  balanced_accuracy = perm_stats
)
write.csv(perm_out,
          "RF_PLSDA_permutation_balanced_accuracy_v4.csv",
          row.names = FALSE)

perm_summary <- data.frame(
  observed_mean_balanced_accuracy = observed_stat,
  permutation_mean = mean(perm_stats, na.rm = TRUE),
  permutation_sd = sd(perm_stats, na.rm = TRUE),
  permutation_p_value = perm_p,
  n_permutations = N_PERM
)
write.csv(perm_summary,
          "RF_PLSDA_permutation_summary_v4.csv",
          row.names = FALSE)

cat("\nPermutation-test summary:\n")
print(perm_summary)

cat("\nLeakage-safe validation complete.\n")

