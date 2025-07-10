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

# Combine 'finger' and 'wrist' into 'hand'
meta_sub$TissueGroup <- as.character(meta_sub$Anatomical_Location)
meta_sub$TissueGroup[meta_sub$TissueGroup %in% c("finger", "wrist")] <- "hand"
meta_sub$TissueGroup <- factor(meta_sub$TissueGroup)

# Select top TFs by variance on **raw** data
top_n <- 250
tf_variances <- apply(data_sub, 1, var)
top_tfs <- names(sort(tf_variances, decreasing = TRUE))[1:top_n]
data_top <- data_sub[top_tfs, ]

# Z-scale across TFs
scaled_top <- t(scale(t(as.matrix(data_top))))
scaled_top <- scaled_top[complete.cases(scaled_top), ]

# Build annotation dataframe
annotation_col <- data.frame(
  Tissue = meta_sub$TissueGroup,
  Treatment = meta_sub$Treatment,
  row.names = samples_to_use
)

# Order columns by tissue
annotation_col$SampleID <- rownames(annotation_col)
annotation_col <- annotation_col %>% arrange(Tissue, Treatment)
ordered_samples <- annotation_col$SampleID
ordered_tissues <- annotation_col$Tissue

# Create gaps between tissue groups
gaps_col <- which(ordered_tissues[-1] != ordered_tissues[-length(ordered_tissues)])

# Reorder scaled data and clean annotation
scaled_top <- scaled_top[, ordered_samples]
rownames(annotation_col) <- annotation_col$SampleID
annotation_col$SampleID <- NULL

# Annotation colors (update for combined 'hand')
ann_colors <- list(
  Tissue = c("hand" = "#E69F00", "hip" = "#56B4E9", "knee" = "#F0E442"),
  Treatment = c("Control" = "#009E73", "medium" = "#009E73", "TNF" = "#D55E00")
)
my_colors <- colorRampPalette(c("navy", "skyblue", "lightgreen"))(100)

# Plot heatmap
pheatmap(scaled_top,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_colnames = FALSE,
         show_rownames = FALSE,
         fontsize = 12,
         color = my_colors,
         gaps_col = gaps_col,
         main = "Top TFs - All Samples")
