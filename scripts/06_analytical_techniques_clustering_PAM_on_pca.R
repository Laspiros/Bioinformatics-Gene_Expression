# ===========================================================================
#  Analytical Techniques - Clustering - PAM (K-MEDOIDS)
# ===========================================================================

# Run the PAM algorithm, asking for 4 groups
# Note: With 4,000+ genes, this might take a few seconds to calculate!
pam_result <- pam(spellman_matrix, k = 4)

# Attach the PAM labels to yhe dataframe
pca_genes_df$PAM_Cluster <- as.factor(pam_result$clustering)

# Plot the PAM clusters on the PCA map
pca_genes_df %>%
  tidyplot(x = PC1, y = PC2, color = PAM_Cluster) %>%
  add_data_points(alpha = 0.8, size = 1.5) %>%
  adjust_size(width = 140, height = 140)

# Plot with limits lines
# xtract the medoids from the PAM result to define the region centers
# and ensure boundaries are correctly centered around the perfect representative.

medoids <- as.data.frame(pam_result$medoids)

# Create the plot with direct ggplot2 and ggforce
ggplot(pca_genes_df, aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(alpha = 0.8, size = 1.5) + # Draws the individual data points
  # Add the ggforce convex hulls with outlines for discrete boundaries
  geom_mark_hull(aes(fill = Cluster, label = Cluster),
                 concavity = 10, # Adjusts how tightly the boundary line fits
                 expand = unit(1, "mm"), # Pushes the line slightly out from the points
                 alpha = 0.2, # Light transparent fill
                 color = "black", # Clear black border line
                 size = 1 # Thickness of the border
  ) +
  scale_fill_discrete(guide = "none") + # Hides the separate legend for the fill
  theme_minimal() + # Use a clean theme
  theme(legend.position = "right") +
  ggtitle("Yeast PCA: PAM Clusters with Discrete Boundary Lines")

# ==============================================================================
# 3-WAY CROSS-TABULATION MATRIX (K-MEANS vs HC vs PAM)
# ==============================================================================

# Ensure the PAM cluster column name matches what is in the dataframe 
three_way_table <- table(
  KMeans = pca_genes_df$Cluster, 
  Hierarchical = pca_genes_df$HC_Cluster,
  PAM = pca_genes_df$PAM_Cluster
)

# Print the 3D table flattened out to view the overlaps
cat("\n--- 3-Way Flattened Intersection Matrix ---\n")
print(ftable(three_way_table))

