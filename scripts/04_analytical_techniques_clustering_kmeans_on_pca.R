# ===========================================================================
#  Analytical Techniques - Clustering - K-Means
# ===========================================================================

# K-MEANS CLUSTERING
set.seed(42) 
kmeans_result <- kmeans(spellman_matrix, centers = 4)

# Attach the new cluster labels directly to PCA dataframe
pca_genes_df$Cluster <- as.factor(kmeans_result$cluster)

# Visualize the Clusters on  PC1 vs PC2 map
pca_genes_df %>%
  tidyplot(x = PC1, y = PC2, color = Cluster) %>%
  add_data_points(alpha = 0.8, size = 1.5) %>%
  adjust_size(width = 140, height = 140)

# Plot the data with boundary lines
ggplot(pca_genes_df, aes(x = PC1, y = PC2, color = Cluster)) +
  # Draw the dots (made slightly transparent so we can see the lines better)
  geom_point(alpha = 0.4, size = 1.5) +
  
  # Draw the boundary lines! 
  # level = 0.95 means the line encloses the core 95% of the cluster
  stat_ellipse(linewidth = 1.5, linetype = "solid") +
  
  # Clean up the visual style
  theme_minimal() +
  labs(
    title = "K-Means Clusters with Boundary Lines",
    subtitle = "95% Confidence Ellipses mapped on PCA Space",
    x = "Principal Component 1",
    y = "Principal Component 2"
  )

# VISUALIZATION: TIME-SERIES RIBBON PLOT
# Create a small table linking each Gene ID to its new K-Means Cluster
clusters_df <- data.frame(
  time = rownames(spellman_matrix),
  Cluster = as.factor(kmeans_result$cluster)
)
# Merge these cluster labels into  'spellman_long' timeline dataset
spellman_clustered <- spellman_long %>%
  left_join(clusters_df, by = "time")

# Generate the Ribbon Line Plot using tidyplots
spellman_clustered %>%
  tidyplot(x = time_point, y = expression, color = Cluster) %>%
  add_mean_line(linewidth = 1.5) %>%    # Draws the solid trend line for the cluster average
  add_sem_ribbon(alpha = 0.3) %>%       # Draws the shaded ribbon (variance) around the line
  split_plot(by = Cluster) %>%          # Splits them into 4 separate mini-plots for extreme clarity!
  adjust_size(width = 200, height = 80)

# Gene Expression Smoothing
# Sort the data chronologically and calculate a Rolling Average
# We use a 'window' of 3 time points. It averages X40, X50, and X60 together 
# to calculate the true midpoint, destroying the random noise!
spellman_smoothed <- spellman_clustered %>%
  group_by(Cluster) %>%
  arrange(time_point) %>%
  mutate(
    smoothed_expression = rollmean(expression, k = 3, fill = NA, align = "center")
  ) %>%
  ungroup()

# Plot the final, smoothed trajectories!
spellman_smoothed %>%
  filter(!is.na(smoothed_expression)) %>% # Removes the empty edges caused by the rolling window
  tidyplot(x = time_point, y = smoothed_expression, color = Cluster) %>%
  add_line(linewidth = 2) %>%           # Draws the beautifully smoothed line
  split_plot(by = Cluster) %>%
  adjust_size(width = 200, height = 80)

# Split and save genes ID's for each cluster
# Extract Gene IDs from row names and split them by Cluster label
cluster_gene_lists <- split(rownames(pca_genes_df), pca_genes_df$Cluster)

# Assign them to distinct standalone vector groups in R workspace
cluster1_genes <- cluster_gene_lists$`1`
cluster2_genes <- cluster_gene_lists$`2`
cluster3_genes <- cluster_gene_lists$`3`
cluster4_genes <- cluster_gene_lists$`4`

# Quick print verification to see how many genes landed in each cluster
cat("--- Gene Counts per Cluster ---\n")
cat("Cluster 1: ", length(cluster1_genes), " genes\n", sep = "")
cat("Cluster 2: ", length(cluster2_genes), " genes\n", sep = "")
cat("Cluster 3: ", length(cluster3_genes), " genes\n", sep = "")
cat("Cluster 4: ", length(cluster4_genes), " genes\n", sep = "")

# Export the groups directly to separate CSV files
write.csv(data.frame(Gene_ID = cluster1_genes), "cluster1_genes.csv", row.names = FALSE)
write.csv(data.frame(Gene_ID = cluster2_genes), "cluster2_genes.csv", row.names = FALSE)
write.csv(data.frame(Gene_ID = cluster3_genes), "cluster3_genes.csv", row.names = FALSE)
write.csv(data.frame(Gene_ID = cluster4_genes), "cluster4_genes.csv", row.names = FALSE)

cat("\nSuccess: All four cluster groups have been saved as CSV files!\n")


