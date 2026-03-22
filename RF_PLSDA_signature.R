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

############################################################
## RF + PLS-DA + Venn + Heatmap of discriminant signature ##
############################################################

## 0. Packages ----
pkgs <- c(
  "dplyr", "tibble", "matrixStats",
  "randomForest", "mixOmics",
  "VennDiagram", "pheatmap"
)
missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(missing)) {
  install.packages(setdiff(missing, c("mixOmics","randomForest","VennDiagram")))
  if ("mixOmics" %in% missing || "randomForest" %in% missing || "VennDiagram" %in% missing) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(setdiff(missing, pkgs[!pkgs %in% c("mixOmics","randomForest","VennDiagram")]), ask = FALSE)
  }
}
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
expr_mat_log <- log1p(expr_mat)

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
grDevices::dev.copy(pdf, file = "RF_PLSDA_venn.pdf", width = 10, height = 10)
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
                     filename = "signature_heatmap.pdf",  # ← Uses correct device + handles dev.off()
                     width = 8,   # inches
                     height = 10  # inches
                     )
  
  
  
  
} else {
  cat("Signature gene set is too small for a meaningful heatmap.\n")
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
 write.csv(signature_table, "signature_genes_RF_PLSDA.csv", row.names = FALSE)

cat("Pipeline complete.\n")

