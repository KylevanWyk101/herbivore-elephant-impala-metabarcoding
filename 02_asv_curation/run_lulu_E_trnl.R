# Load required packages
if (!requireNamespace("lulu", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("vegan", quietly = TRUE)) install.packages("vegan")
devtools::install_github("tobiasgf/lulu")

library(lulu)

# Set working directory (adjust if needed)
setwd("C:/University/2025/kyle_research_project/Data_processing/Seq_data/Manifest_qiime2_files/E_trnl_output/LULU_E_trnl")

readLines("E_trnl_feature_table_min50.tsv", n = 5)
# Step 1: Load the ASV feature table
# go into the OTU tsv file and manually remove the '#' TO ENSURE THE 2ND LINE IS THE HEADER OF IMPORTED TABLE
asv_table <- read.table("E_trnl_feature_table_min50.tsv",
												sep = "\t",
												header = TRUE,
												row.names = 1,
												comment.char = "#",
												check.names = FALSE)


# Step 2: Transpose the table (samples must be columns, ASVs rows)
# asv_table <- t(asv_table)

# Step 3: Load the match list file
match_list <- read.table("trnl_E_lulu_match_list.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)

# Step 4: Run LULU
lulu_result <- lulu(asv_table, match_list)

# Step 5: Save curated ASV table
write.table(lulu_result$curated_table,
						file = "E_trnl_LULU_curated_table.tsv",
						sep = "\t", quote = FALSE)

# Step 6: Save ASV ID tracking (to see what was merged or kept)
write.table(lulu_result$curated_ids,
						file = "E_trnl_LULU_curation_map.tsv",
						sep = "\t", quote = FALSE, row.names = FALSE)

# Optional: Print summary
cat("Initial number of ASVs:", nrow(asv_table), "\n")
cat("Number of curated ASVs after LULU:", nrow(lulu_result$curated_table), "\n")
