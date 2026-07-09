# ==============================================================================
# EXTRACTING SUPER-CORE CONSENSUS GENES
# ==============================================================================

# Extract PAM lists just like before
pam_gene_lists <- split(rownames(pca_genes_df), pca_genes_df$PAM_Cluster)
pam_cluster1_genes <- pam_gene_lists$`1`
# (Repeat splitting for pam_cluster2, 3, 4 vectors as needed)

# Example: Finding genes that ALL THREE models agree belong to the primary groups
# (Adjust cluster numbers based on your specific matrix alignment)
super_core_genes <- intersect(
  intersect(cluster2_genes, hc_cluster2_genes), # KMeans 2 & HC 2
  pam_cluster1_genes                             # Cross-reference with matching PAM cluster
)

cat("\nIsolated", length(super_core_genes), "Super-Core Consensus Genes approved by all 3 algorithms.\n")
write.csv(data.frame(Gene_ID = super_core_genes), "super_core_consensus.csv", row.names = FALSE)

# ==============================================================================
# AUTOMATED 3-WAY SUPER-CORE EXTRACTION (ERROR-PROOF)
# ==============================================================================

# Ensure the list splits are freshly active
cluster_gene_lists <- split(rownames(pca_genes_df), pca_genes_df$Cluster)
hc_gene_lists      <- split(rownames(pca_genes_df), pca_genes_df$HC_Cluster)
pam_gene_lists     <- split(rownames(pca_genes_df), pca_genes_df$PAM_Cluster)

# Extract Super-Core Phase 1 (From Row 19: 793 genes)
phase_1_core <- intersect(intersect(cluster_gene_lists$`3`, hc_gene_lists$`1`), pam_gene_lists$`2`)

# Extract Super-Core Phase 2 (From Row 38: 739 genes)
phase_2_core <- intersect(intersect(cluster_gene_lists$`2`, hc_gene_lists$`2`), pam_gene_lists$`3`)

# Extract Super-Core Phase 3 (From Row 49: 628 genes)
phase_3_core <- intersect(intersect(cluster_gene_lists$`1`, hc_gene_lists$`1`), pam_gene_lists$`4`)

# Extract Super-Core Phase 4 (From Row 53: 540 genes)
phase_4_core <- intersect(intersect(cluster_gene_lists$`1`, hc_gene_lists$`2`), pam_gene_lists$`4`)

# Extract Super-Core Phase 5 (From Row 6: 253 genes)
phase_5_core <- intersect(intersect(cluster_gene_lists$`2`, hc_gene_lists$`2`), pam_gene_lists$`1`)


# --- Print Total Consensus Mapping Verification ---
cat("\n--- Corrected Super-Core Genome Mapping ---\n")
cat("Super-Core Phase 1 Matrix Yield: ", length(phase_1_core), " genes\n", sep = "")
cat("Super-Core Phase 2 Matrix Yield: ", length(phase_2_core), " genes\n", sep = "")
cat("Super-Core Phase 3 Matrix Yield: ", length(phase_3_core), " genes\n", sep = "")
cat("Super-Core Phase 4 Matrix Yield: ", length(phase_4_core), " genes\n", sep = "")
cat("Super-Core Phase 5 Matrix Yield: ", length(phase_5_core), " genes\n", sep = "")


# --- Save all verified elite subsets safely to your working directory ---
write.csv(data.frame(Gene_ID = phase_1_core), "super_core_phase_1.csv", row.names = FALSE)
write.csv(data.frame(Gene_ID = phase_2_core), "super_core_phase_2.csv", row.names = FALSE)
write.csv(data.frame(Gene_ID = phase_3_core), "super_core_phase_3.csv", row.names = FALSE)
write.csv(data.frame(Gene_ID = phase_4_core), "super_core_phase_4.csv", row.names = FALSE)
write.csv(data.frame(Gene_ID = phase_5_core), "super_core_phase_5.csv", row.names = FALSE)

cat("\nSuccess: All high-yield super-core files saved to your project directory!\n")
