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

# Filter for Control and Medium
meta_ctrl <- meta_sub[meta_sub$Treatment %in% c("Control", "medium"), ]
ctrl_samples <- meta_ctrl$Sample.ID
data_ctrl <- data_sub[, ctrl_samples]

# Combine finger and wrist → hand
meta_ctrl$TissueGroup <- as.character(meta_ctrl$Anatomical_Location)
meta_ctrl$TissueGroup[meta_ctrl$TissueGroup %in% c("finger", "wrist")] <- "hand"
meta_ctrl$TissueGroup <- factor(meta_ctrl$TissueGroup)

# Top TFs by variance on raw data
top_n <- 250
tf_variances <- apply(data_sub, 1, var)
top_tfs <- names(sort(tf_variances, decreasing = TRUE))[1:top_n]
data_top <- data_sub[top_tfs, ]
scaled_top <- t(scale(t(as.matrix(data_top))))
scaled_top <- scaled_top[complete.cases(scaled_top), ]
scaled_ctrl <- scaled_top[, ctrl_samples]

# Annotations and ordering
annotation_col <- data.frame(
  Tissue = meta_ctrl$TissueGroup,
  Treatment = meta_ctrl$Treatment,
  row.names = ctrl_samples
)
annotation_col$SampleID <- rownames(annotation_col)
annotation_col <- annotation_col %>% arrange(Tissue, Treatment)
ordered_samples <- annotation_col$SampleID
ordered_tissues <- annotation_col$Tissue
gaps_col <- which(ordered_tissues[-1] != ordered_tissues[-length(ordered_tissues)])
scaled_ctrl <- scaled_ctrl[, ordered_samples]
rownames(annotation_col) <- annotation_col$SampleID
annotation_col$SampleID <- NULL

# Colors
ann_colors <- list(
  Tissue = c("hand" = "#E69F00", "hip" = "#56B4E9", "knee" = "#F0E442"),
  Treatment = c("Control" = "#009E73", "medium" = "#009E73")
)
my_colors <- colorRampPalette(c("navy", "skyblue", "lightgreen"))(100)

# Plot
pheatmap(scaled_ctrl,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_colnames = FALSE,
         show_rownames = FALSE,
         fontsize = 12,
         color = my_colors,
         gaps_col = gaps_col,
         main = "Top TFs - Control/Medium")
