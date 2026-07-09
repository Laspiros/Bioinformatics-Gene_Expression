# ==============================================================================
# GMM SOFT CLUSTERING ON UMAP TOPOLOGY
# ==============================================================================

library(mclust)
library(ggplot2)

cat("Running GMM on the UMAP projection...\n")

# Isolate the 2D UMAP coordinates we just generated
umap_pca_coords <- pca_genes_df[, c("UMAP_PCA_1", "UMAP_PCA_2")]

# Run the Gaussian Mixture Model (k=4 phases)
set.seed(42)
gmm_umap_model <- Mclust(umap_pca_coords, G = 4)

# Attach the primary assignments and the UNCERTAINTY scores to the dataframe
pca_genes_df$GMM_UMAP_Cluster <- as.factor(gmm_umap_model$classification)
pca_genes_df$GMM_UMAP_Uncertainty <- gmm_umap_model$uncertainty

# ==============================================================================
# VISUALIZATION: THE UNCERTAINTY MAP
# ==============================================================================

# Create a custom label to highlight genes with highly split probabilities (Uncertainty > 0.3)
pca_genes_df$Is_Bridge <- ifelse(pca_genes_df$GMM_UMAP_Uncertainty > 0.3, "Transient Bridge", "Core Gene")

# Plot the UMAP Manifold, coloring the cores by phase, and highlighting the bridges in Red
ggplot(pca_genes_df, aes(x = UMAP_PCA_1, y = UMAP_PCA_2)) +
  # Plot the stable core genes with soft transparency
  geom_point(data = subset(pca_genes_df, Is_Bridge == "Core Gene"), 
             aes(color = GMM_UMAP_Cluster), alpha = 0.4, size = 1.2) +
  # Overlay the highly uncertain "Bridge Genes" in stark red
  geom_point(data = subset(pca_genes_df, Is_Bridge == "Transient Bridge"), 
             color = "red", alpha = 0.9, size = 2) +
  theme_minimal() +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "GMM Probabilistic Mapping on UMAP Manifold",
    subtitle = "Red points indicate highly transient bridge genes caught between phase boundaries",
    x = "UMAP Dimension 1",
    y = "UMAP Dimension 2",
    color = "Phase Assignment"
  )

# ==============================================================================
# EXTRACT THE UMAP-GMM BRIDGE GENES
# ==============================================================================

umap_bridge_genes <- rownames(pca_genes_df[pca_genes_df$Is_Bridge == "Transient Bridge", ])
cat("\nIsolated", length(umap_bridge_genes), "transient bridge genes across UMAP boundaries.\n")

write.csv(data.frame(Gene_ID = umap_bridge_genes), "umap_gmm_transient_bridges.csv", row.names = FALSE)
