# ===========================================================================
#  Analytical Techniques - Clustering - Hierarchical Clustering
# ===========================================================================

# Calculate the actual mathematical distance between every single gene
distance_matrix <- dist(spellman_matrix, method = "euclidean")

# Build the family tree using Ward's Method
hc_result <- hclust(distance_matrix, method = "ward.D2")

# "Cut" the tree to create 4 clusters to compare with your K-Means result
pca_genes_df$HC_Cluster <- as.factor(cutree(hc_result, k = 4))

# Plot the new Hierarchical Clusters on the PCA map
pca_genes_df %>%
  tidyplot(x = PC1, y = PC2, color = HC_Cluster) %>%
  add_data_points(alpha = 0.8, size = 1.5) %>%
  adjust_size(width = 180, height = 180)

# Plot the Hierarchical data with boundary lines
ggplot(pca_genes_df, aes(x = PC1, y = PC2, color = HC_Cluster)) +
  # Draw the dots (made transparent to highlight the lines)
  geom_point(alpha = 0.4, size = 1.5) +
  
  # Draw the boundary lines (enclosing the core 95% of each hierarchical family)
  stat_ellipse(linewidth = 1.5, linetype = "solid") +
  
  # Clean up the visual style
  theme_minimal() +
  labs(
    title = "Hierarchical Clustering with Boundary Lines",
    subtitle = "Ward's Method Branches mapped on PCA Space",
    x = "Principal Component 1",
    y = "Principal Component 2",
    color = "HC Branch"
  )

# SPLIT AND SAVE GENE IDs FOR HIERARCHICAL CLUSTERING (HC)

# Extract Gene IDs from row names and split them by the HC Cluster labels
hc_gene_lists <- split(rownames(pca_genes_df), pca_genes_df$HC_Cluster)

# Assign them to distinct standalone vectors 
hc_cluster1_genes <- hc_gene_lists$`1`
hc_cluster2_genes <- hc_gene_lists$`2`
hc_cluster3_genes <- hc_gene_lists$`3`
hc_cluster4_genes <- hc_gene_lists$`4`

# Print counts to verify distribution
cat("--- HC Gene Counts per Branch ---\n")
cat("HC Cluster 1: ", length(hc_cluster1_genes), " genes\n", sep = "")
cat("HC Cluster 2: ", length(hc_cluster2_genes), " genes\n", sep = "")
cat("HC Cluster 3: ", length(hc_cluster3_genes), " genes\n", sep = "")
cat("HC Cluster 4: ", length(hc_cluster4_genes), " genes\n", sep = "")

# Export the HC groups directly to separate CSV files
write.csv(data.frame(Gene_ID = hc_cluster1_genes), "hc_cluster1_genes.csv", row.names = FALSE)
write.csv(data.frame(Gene_ID = hc_cluster2_genes), "hc_cluster2_genes.csv", row.names = FALSE)
write.csv(data.frame(Gene_ID = hc_cluster3_genes), "hc_cluster3_genes.csv", row.names = FALSE)
write.csv(data.frame(Gene_ID = hc_cluster4_genes), "hc_cluster4_genes.csv", row.names = FALSE)


# ALGORITHM COMPARISON & INTERSECTION MATRIX
# Create a cross-tabulation table comparing the two clustering vectors
comparison_matrix <- table(KMeans = pca_genes_df$Cluster, Hierarchical = pca_genes_df$HC_Cluster)

# Print the matrix directly to console
cat("\n--- Intersection Matrix (K-Means vs Hierarchical) ---\n")
print(comparison_matrix)

# ==============================================================================
# EXTRACTING TRUE ENSEMBLE CONSENSUS GROUPS (BASED ON RESULTS matrix)
# ==============================================================================

# Consensus Group A: The Massive Cluster 2 / Branch 2 Alignment (1,042 genes)
consensus_core_A <- intersect(cluster2_genes, hc_cluster2_genes)

# Consensus Group B: The Massive Cluster 3 / Branch 1 Alignment (1,041 genes)
consensus_core_B <- intersect(cluster3_genes, hc_cluster3_genes)

# Consensus Group C: The Cluster 1 / Branch 1 Split Vector (794 genes)
consensus_split_C1_HC1 <- intersect(cluster1_genes, hc_cluster1_genes)

# Consensus Group D: The Cluster 1 / Branch 2 Split Vector (786 genes)
consensus_split_C1_HC2 <- intersect(cluster1_genes, hc_cluster2_genes)


# --- Print Verification Reports ---
cat("\n--- Extracted Consensus Gene Counts ---\n")
cat("Core Consensus Group A (K2 & HC2): ", length(consensus_core_A), " genes\n", sep = "")
cat("Core Consensus Group B (K3 & HC1): ", length(consensus_core_B), " genes\n", sep = "")
cat("Split Boundary Group C (K1 & HC1): ", length(consensus_split_C1_HC1), " genes\n", sep = "")
cat("Split Boundary Group D (K1 & HC2): ", length(consensus_split_C1_HC2), " genes\n", sep = "")


# --- Export Data to Workspace Directory ---
write.csv(data.frame(Gene_ID = consensus_core_A), "core_consensus_group_A.csv", row.names = FALSE)
write.csv(data.frame(Gene_ID = consensus_core_B), "core_consensus_group_B.csv", row.names = FALSE)
write.csv(data.frame(Gene_ID = consensus_split_C1_HC1), "boundary_consensus_group_C.csv", row.names = FALSE)
write.csv(data.frame(Gene_ID = consensus_split_C1_HC2), "boundary_consensus_group_D.csv", row.names = FALSE)

cat("\nSuccess: Clean consensus groups exported for downstream biological validation!\n")
