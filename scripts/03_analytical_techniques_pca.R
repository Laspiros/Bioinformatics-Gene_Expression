# ===========================================================================
#  Analytical Techniques - PCA 
# ===========================================================================

# Compute PCA on the genes (rows are observations, columns are features)
# scale. = TRUE standardizes the time points so they carry equal weight
pca_result <- prcomp(spellman_matrix, scale. = TRUE)

# Extract the coordinate scores for the genes (found in the 'x' slot)
pca_genes_df <- as.data.frame(pca_result$x)

# Create the PCA Scatter Plot using tidyplots
pca_genes_df %>%
  tidyplot(x = PC1, y = PC2) %>%
  add_data_points() %>%
  adjust_size(width = 120, height = 120)

# PCA Scree plot
# Building and Plotting the Scree Data
# Calculate the percentage of variance explained by each principal component
pca_var <- pca_result$sdev^2
pca_var_perc <- round(pca_var / sum(pca_var) * 100, 1)

# Create the 'scree_data' dataframe (looking at just the top 10 PCs)
scree_data <- data.frame(
  PC = factor(paste0("PC", 1:10), levels = paste0("PC", 1:10)), 
  Variance = pca_var_perc[1:10]
)

# Generate the Scree Plot using tidyplots
scree_data %>%
  tidyplot(x = PC, y = Variance) %>%
  add_line(linewidth = 1, group = 1) %>%  # group = 1 ensures the line connects categorical X-axis points
  add_data_points(size = 3) %>%           # Adds dots at the exact values
  adjust_size(width = 160, height = 100)

# Mapping PC3 to Color Gradient
pca_genes_df %>%
  tidyplot(x = PC1, y = PC2, color = PC3) %>%
  add_data_points(alpha = 0.8, size = 1.5) %>% # alpha adds slight transparency
  adjust_size(width = 140, height = 140)

# Plotting PC1 directly against PC3
pca_genes_df %>%
  tidyplot(x = PC1, y = PC3) %>%
  add_data_points(alpha = 0.5, size = 1.5) %>% 
  adjust_size(width = 120, height = 120)

