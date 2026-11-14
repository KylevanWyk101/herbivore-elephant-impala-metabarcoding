# -------------------------------
# Differential Abundance: trnl dataset (Genus-level, DESeq2)
# -------------------------------

# Load required libraries
library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggrepel)

# Set working directory
setwd("C:/University/2025/kyle_research_project/Statistics/Differential_abundance/trnl/DESeq2_DA_test")

# Read the combined dataset
data <- read.csv("trnl_combined_dietlist.csv", stringsAsFactors = FALSE)

# Rename 'Life form' to 'Life_form' to avoid issues with spaces
colnames(data)[colnames(data) == "Life form"] <- "Life_form"

# Define taxonomy columns
taxonomy_cols <- c("ASV_ID", "Family", "Genus", "Species", "Life_form")

# Identify sample columns (exclude taxonomy columns)
sample_cols <- setdiff(colnames(data), taxonomy_cols)

# Verify that sample columns are numeric
non_numeric_cols <- sample_cols[!sapply(data[sample_cols], is.numeric)]
if (length(non_numeric_cols) > 0) {
	cat("Warning: The following sample columns are non-numeric and will be excluded:\n")
	cat(paste(non_numeric_cols, collapse = ", "), "\n")
	sample_cols <- setdiff(sample_cols, non_numeric_cols)
}

# Collapse data by Genus (sum sample counts for each genus)
genus_data <- data %>%
	group_by(Genus) %>%
	summarise(across(all_of(sample_cols), sum, na.rm = TRUE)) %>%
	as.data.frame()

# Set Genus as row names for DESeq2
rownames(genus_data) <- genus_data$Genus
genus_data <- genus_data[, sample_cols, drop = FALSE]


# Report dataset dimensions
cat("Number of samples:", ncol(genus_data), "\n")
cat("Elephant samples:", sum(grepl("_E$", colnames(genus_data))), "\n")
cat("Impala samples:", sum(grepl("_I$", colnames(genus_data))), "\n")

# Create sample metadata
sample_names <- colnames(genus_data)
condition <- ifelse(grepl("_E$", sample_names), "Elephant", "Impala")
col_data <- data.frame(
	sample = sample_names,
	condition = factor(condition, levels = c("Impala", "Elephant")),
	row.names = sample_names
)

# Create DESeq2
dds <- DESeqDataSetFromMatrix(
	countData = genus_data,
	colData = col_data,
	design = ~ condition
)

# Use poscounts for size factor estimation to handle sparse data
dds <- estimateSizeFactors(dds, type = "poscounts")

# Run DESeq2 analysis
dds <- DESeq(dds)

# Extract results (Elephant vs. Impala)
res <- results(dds, contrast = c("condition", "Elephant", "Impala"))
res <- as.data.frame(res)

# Add Genus column and handle NA p-values
res$Genus <- rownames(res)
res$padj <- ifelse(is.na(res$padj), 1, res$padj)  # Replace NA padj with 1

# Mark significant genera (THIS WAS MISSING IN SOME RUNS)
res$significant <- res$padj < 0.05 & abs(res$log2FoldChange) > 1

# Save DESeq2 results to CSV
write.csv(res, "deseq_results.csv", row.names = FALSE)
cat("DESeq2 results saved as 'deseq_results.csv'\n")

# Volcano plot with ggrepel to prevent overlaps and keep labels above axis
volcano_plot <- ggplot(res, aes(x = log2FoldChange, y = -log10(padj))) +
	geom_point(aes(color = significant), size = 2.5) +
	scale_color_manual(values = c("grey50", "red"),
										 labels = c("Not Significant", "Significant")) +
	geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
	geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
	geom_text_repel(
		data = subset(res, significant),
		aes(label = Genus),
		size = 4,
		color = "black",
		max.overlaps = 100,              
		box.padding = 1.0,              # Increased from 0.6 for more space around labels
		point.padding = 0.5,            # Increased from 0.3 for more space around points
		segment.color = "grey60",
		segment.size = 0.2,             # Reduced line thickness
		ylim = c(0, 80),                # keep labels within y range
		direction = "both",             # Allow both x and y direction for better placement
		force = 3,                      # to push labels apart more aggressively
		seed = 42                       # Set seed for reproducibility
	) +
	coord_cartesian(ylim = c(0, 80), clip = "off") +  # don’t cut off labels
	annotate("text", x = -35, y = 75, label = "Impala (Higher Abundance)",
					 color = "black", size = 5, hjust = 0, fontface = "bold") +
	annotate("text", x = 35, y = 75, label = "Elephant (Higher Abundance)",
					 color = "black", size = 5, hjust = 1, fontface = "bold") +
	theme_minimal(base_size = 14) +
	theme(
		panel.background = element_rect(fill = "white", color = NA),
		plot.background = element_rect(fill = "white", color = NA),
		panel.grid.major = element_line(color = "grey85", size = 0.3),
		panel.grid.minor = element_line(color = "grey92", size = 0.2),
		axis.text = element_text(color = "black", size = 13),
		axis.title = element_text(color = "black", size = 13),
		legend.text = element_text(color = "black", size = 13),
		legend.title = element_text(color = "black", size = 14, face = "bold"),
		plot.title = element_text(color = "black", size = 18, face = "bold", hjust = 0.5),
		plot.margin = margin(20, 20, 60, 20)  # Increased top margin to 60
	) +
	labs(
		title = "Volcano Plot: trnl Elephant vs. Impala (Genus-Level)",
		x = "Log2 Fold Change (Elephant/Impala)",
		y = expression(-log[10]~"Adjusted P-value"),
		color = "Significance"
	)

# Save volcano plot with increased width
ggsave("trnl_volcano_plot.png", volcano_plot,
			 width = 14, height = 8, dpi = 300, bg = "white")
cat("Volcano plot saved as 'trnl_volcano_plot.png'\n")


