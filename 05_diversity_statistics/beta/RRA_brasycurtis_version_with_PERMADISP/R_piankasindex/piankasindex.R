# =====================================================
# R script: Compute Pianka's Index with pgirmess
# =====================================================
setwd("C:/University/2025/kyle_research_project/Statistics/Beta_diversity/RRA_brasycurtis_version/R_piankasindex")

library(pgirmess)
library(dplyr)

datasets <- list(
	rbcL = c("E_rbcl_diet_rra_only.csv", "I_rbcl_diet_rra_only.csv"),
	trnL = c("E_trnl_diet_rra_only.csv", "I_trnl_diet_rra_only.csv")
)

for (locus in names(datasets)) {
	cat("\n===== Processing locus:", locus, "=====\n")
	files <- datasets[[locus]]
	
	# ---------- Elephant ----------
	df1 <- read.csv(files[1], stringsAsFactors = FALSE)
	rownames(df1) <- df1$ASV_ID
	rra_cols1 <- grep("_RRA$", colnames(df1), value = TRUE)
	abund1 <- t(df1[, rra_cols1, drop = FALSE])
	rownames(abund1) <- gsub("_RRA$", "", rra_cols1)   # sample names
	# colnames(abund1) = ASV_IDs automatically
	
	# ---------- Impala ----------
	df2 <- read.csv(files[2], stringsAsFactors = FALSE)
	rownames(df2) <- df2$ASV_ID
	rra_cols2 <- grep("_RRA$", colnames(df2), value = TRUE)
	abund2 <- t(df2[, rra_cols2, drop = FALSE])
	rownames(abund2) <- gsub("_RRA$", "", rra_cols2)
	
	# ---------- Align ASV columns ----------
	all_asvs <- union(colnames(abund1), colnames(abund2))
	
	abund1_full <- matrix(0, nrow = nrow(abund1), ncol = length(all_asvs),
												dimnames = list(rownames(abund1), all_asvs))
	abund1_full[, colnames(abund1)] <- abund1
	
	abund2_full <- matrix(0, nrow = nrow(abund2), ncol = length(all_asvs),
												dimnames = list(rownames(abund2), all_asvs))
	abund2_full[, colnames(abund2)] <- abund2
	
	# ---------- Combine ----------
	combined <- rbind(abund1_full, abund2_full)
	
	# ---------- Pianka similarity ----------
	pianka_mat <- piankabio(combined, combined)
	
	write.csv(pianka_mat,
						paste0("pianka_matrix_", locus, "_rra.csv"),
						quote = FALSE)
	
	# ---------- Group summaries ----------
	groups <- c(rep("Elephant", nrow(abund1_full)),
							rep("Impala", nrow(abund2_full)))
	names(groups) <- rownames(combined)
	
	sample_names <- rownames(combined)
	results <- data.frame(Group1=character(),
												Group2=character(),
												Pianka=numeric(),
												stringsAsFactors=FALSE)
	
	for (i in 1:(length(sample_names)-1)) {
		for (j in (i+1):length(sample_names)) {
			s1 <- sample_names[i]; s2 <- sample_names[j]
			g1 <- groups[s1]; g2 <- groups[s2]
			val <- pianka_mat[s1, s2]
			results <- rbind(results,
											 data.frame(Group1=g1, Group2=g2, Pianka=val,
											 					 stringsAsFactors=FALSE))
		}
	}
	
	summary <- results %>%
		group_by(Group1, Group2) %>%
		summarise(Pianka = mean(Pianka), .groups="drop")
	
	write.csv(summary,
						paste0("pianka_summary_", locus, "_rra.csv"),
						row.names=FALSE, quote=FALSE)
	
	cat("Saved Pianka outputs for", locus, "\n")
}

cat("\nAll loci processed successfully.\n")

