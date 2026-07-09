# ===========================================================================
#  Analytical Techniques - Statistical testing
# ===========================================================================

# Statistical testing (ANOVA)
# Filter the dataset to look at ONE highly active time point (e.g., minute 130)
spellman_clustered %>%
  filter(time_point == "X130") %>%
  
  # Set up the plot comparing the 4 Clusters
  tidyplot(x = Cluster, y = expression, color = Cluster) %>%
  
  # Add boxplots to show the distribution spread
  add_boxplot(alpha = 0.5) %>%
  
  # Automatically run statistical testing and add the p-value!
  add_test_pvalue() %>% 
  
  adjust_size(width = 160, height = 120)

# Extracting the Gene Names

# Extract the gene IDs for a specific cluster (e.g., PAM Cluster 2)
cluster_2_genes <- pca_genes_df %>%
  rownames_to_column(var = "Gene_ID") %>%   # Converts the hidden row names into a real column
  filter(PAM_Cluster == "2") %>%            # Keeps only the genes in Cluster 2
  pull(Gene_ID)                             # Extracts them as a clean list

# Print the first 30 genes to the console
head(cluster_2_genes, 30)

# Save the full list to a CSV file to share 
write.csv(cluster_2_genes, "Cluster_2_Genes.csv", row.names = FALSE)
