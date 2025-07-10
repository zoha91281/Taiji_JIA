library(pheatmap)
library(readr)
library(dplyr)

# Load data
data <- read.csv("~/Downloads/raw.pageranks.csv", row.names = 1, check.names = FALSE)
meta <- read.csv("~/Downloads/META.CSV", row.names = 1, check.names = FALSE)

# Clean and match samples
meta$Sample.ID <- trimws(meta$Sample.ID)
samples_to_use <- intersect(meta$Sample.ID, colnames(data))
data_sub <- data[, samples_to_use]
meta_sub <- meta[match(samples_to_use, meta$Sample.ID), ]

# Filter for TNF only
meta_tnf <- meta_sub[meta_sub$Treatment == "TNF", ]
tnf_samples <- meta_tnf$Sample.ID
data_tnf <- data_sub[, tnf_samples]

# Combine finger and wrist → hand
meta_tnf$TissueGroup <- as.character(meta_tnf$Anatomical_Location)
meta_tnf$TissueGroup[meta_tnf$TissueGroup %in% c("finger", "wrist")] <- "hand"
meta_tnf$TissueGroup <- factor(meta_tnf$TissueGroup)

# Top TFs by variance on raw data
top_n <- 250
tf_variances <- apply(data_sub, 1, var)
top_tfs <- names(sort(tf_variances, decreasing = TRUE))[1:top_n]
data_top <- data_sub[top_tfs, ]
scaled_top <- t(scale(t(as.matrix(data_top))))
scaled_top <- scaled_top[complete.cases(scaled_top), ]
scaled_tnf <- scaled_top[, tnf_samples]

# Annotations and ordering
annotation_col <- data.frame(
  Tissue = meta_tnf$TissueGroup,
  Treatment = meta_tnf$Treatment,
  row.names = tnf_samples
)
annotation_col$SampleID <- rownames(annotation_col)
annotation_col <- annotation_col %>% arrange(Tissue, Treatment)
ordered_samples <- annotation_col$SampleID
ordered_tissues <- annotation_col$Tissue
gaps_col <- which(ordered_tissues[-1] != ordered_tissues[-length(ordered_tissues)])
scaled_tnf <- scaled_tnf[, ordered_samples]
rownames(annotation_col) <- annotation_col$SampleID
annotation_col$SampleID <- NULL

# Colors
ann_colors <- list(
  Tissue = c("hand" = "#E69F00", "hip" = "#56B4E9", "knee" = "#F0E442"),
  Treatment = c("TNF" = "#D55E00")
)
my_colors <- colorRampPalette(c("navy", "skyblue", "lightgreen"))(100)

# Plot
pheatmap(scaled_tnf,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_colnames = FALSE,
         show_rownames = FALSE,
         fontsize = 12,
         color = my_colors,
         gaps_col = gaps_col,
         main = "Top TFs - TNF Only")
