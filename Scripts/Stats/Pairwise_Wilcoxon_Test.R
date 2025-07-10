# Load data
data <- read.csv("~/Downloads/raw.pageranks.csv", row.names = 1, check.names = FALSE)
meta <- read.csv("~/Downloads/META.CSV", row.names = 1, check.names = FALSE)

# Clean sample IDs
meta$Sample.ID <- trimws(meta$Sample.ID)

# Match samples
samples_to_use <- intersect(meta$Sample.ID, colnames(data))
data <- data[, samples_to_use]
meta <- meta[match(samples_to_use, meta$Sample.ID), ]

# Subset Control samples only
meta_ctrl <- meta[meta$Treatment == "Control", ]
ctrl_samples <- meta_ctrl$Sample.ID
data_ctrl <- data[, ctrl_samples]

# Combine finger and wrist into "hand"
meta_ctrl$TissueGroup <- as.character(meta_ctrl$Anatomical_Location)
meta_ctrl$TissueGroup[meta_ctrl$TissueGroup %in% c("finger", "wrist")] <- "hand"
meta_ctrl$TissueGroup <- factor(meta_ctrl$TissueGroup)

# Update sample-to-tissue map
tissue_ctrl <- meta_ctrl$TissueGroup
names(tissue_ctrl) <- meta_ctrl$Sample.ID

# Define comparisons between tissue groups
comparisons <- list(
  c("hand", "hip"),
  c("hand", "knee"),
  c("hip", "knee")
)

# Initialize result list
wilcox_results <- list()

# Loop over TFs and perform pairwise Wilcoxon tests
for (tf in rownames(data_ctrl)) {
  tf_values <- data_ctrl[tf, ]
  tf_pvals <- list()
  
  for (pair in comparisons) {
    group1 <- pair[1]
    group2 <- pair[2]
    
    vals1 <- as.numeric(tf_values[names(tissue_ctrl)[tissue_ctrl == group1]])
    vals2 <- as.numeric(tf_values[names(tissue_ctrl)[tissue_ctrl == group2]])
    
    if (length(vals1) >= 2 && length(vals2) >= 2) {
      test_result <- wilcox.test(vals1, vals2, exact = FALSE)
      tf_pvals[[paste(group1, "vs", group2)]] <- test_result$p.value
    } else {
      tf_pvals[[paste(group1, "vs", group2)]] <- NA
    }
  }
  
  wilcox_results[[tf]] <- tf_pvals
}

# Combine results into a dataframe
wilcox_df <- do.call(rbind, lapply(wilcox_results, function(x) unlist(x)))
wilcox_df <- as.data.frame(wilcox_df)
wilcox_df$TF <- rownames(wilcox_df)
rownames(wilcox_df) <- NULL

# Reorder columns for clarity
wilcox_df <- wilcox_df[, c("TF", "hand vs hip", "hand vs knee", "hip vs knee")]

# View and save
head(wilcox_df)
write.csv(wilcox_df, file = "~/Downloads/control_wilcox_results_grouped_hand.csv", row.names = FALSE)
