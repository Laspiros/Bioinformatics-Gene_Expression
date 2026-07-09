# ==============================================================================
# UMAP PROJECTION ON PCA-COMPRESSED DATA (NOISE-FILTERED TOPOLOGY)
# ==============================================================================

library(umap)
library(ggplot2)

cat("Running UMAP on PCA coordinates. This will be faster and cleaner...\n")

# We use the pca_coords we created earlier (e.g., PC1 to PC3 or up to PC10)
pca_coords <- pca_genes_df[, grep("^PC", colnames(pca_genes_df))]

# Run UMAP directly on the PCA coordinates
set.seed(42) # Lock the random seed for reproducibility
umap_pca_results <- umap(pca_coords)

# Extract the new 2D coordinates and add them to the master dataframe
pca_genes_df$UMAP_PCA_1 <- umap_pca_results$layout[, 1]
pca_genes_df$UMAP_PCA_2 <- umap_pca_results$layout[, 2]

# 4. Visualize the Noise-Filtered UMAP Manifold, coloring by the PAM Clusters
ggplot(pca_genes_df, aes(x = UMAP_PCA_1, y = UMAP_PCA_2, color = PAM_Cluster)) +
  geom_point(alpha = 0.6, size = 1.5) +
  theme_minimal() +
  scale_color_brewer(palette = "Set1") + 
  labs(
    title = "Noise-Filtered UMAP Projection (via PCA)",
    subtitle = "Non-linear manifold extraction of the principal component space",
    x = "UMAP Dimension 1",
    y = "UMAP Dimension 2"
  )

