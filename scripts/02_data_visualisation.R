# ===========================================================================
#  Data Visualization 
# ===========================================================================

# Reshape exact dataset from WIDE to LONG format
spellman_long <- spellman_data %>%
  pivot_longer(
    cols = starts_with("X"), 
    names_to = "time_point", 
    values_to = "expression"
  )

# The Overall Distribution Plot
spellman_long %>%
  tidyplot(x = expression) %>%
  add_histogram(bins = 100) %>%  
  add_reference_lines(x = 0) %>%
  adjust_x_axis_title("Overall Distribution Plot") # <-- Re-labels the X-axis

# Boxplots
# Extract the exact sequential order of time points 
time_order <- colnames(spellman_data)[-1]

# Apply the sequential order to the long dataset so the plot reads left-to-right
spellman_long <- spellman_long %>%
  mutate(time_point = factor(time_point, levels = time_order))

# Generate the fully colored Boxplot
spellman_long %>%
  tidyplot(x = time_point, y = expression, color = time_point) %>%
  add_boxplot() %>%
  remove_legend() %>% 
  adjust_size(width = 160, height = 90)

# Violin Plot
# Generate the Violin Plot 
spellman_long %>%
  tidyplot(x = time_point, y = expression, color = time_point) %>%
  add_violin() %>%        # <-- Just change this line!
  remove_legend() %>% 
  adjust_size(width = 160, height = 90)

# Heatmap with tidyplots
# Select the first 40 unique genes from the dataset so the Y-axis is readable
genes_to_plot <- unique(spellman_long$time)[1:40]

# Filter the long dataset and pipe it into tidyplots
spellman_long %>%
  filter(time %in% genes_to_plot) %>%               # Keep only those 40 genes
  tidyplot(x = time_point, y = time, color = expression) %>%
  add_heatmap() %>%                                 # Draws the color grid
  adjust_size(width = 160, height = 160)            # Make it tall enough for the gene names

# Heatmap with Hierarchical Clustering
rownames(spellman_matrix) <- spellman_data$time

# Generate the Clustered Heatmap
pheatmap(
  mat = spellman_matrix, 
  cluster_rows = TRUE,       # YES: Group similar genes together
  cluster_cols = FALSE,      # NO: Keep time points in chronological order
  show_rownames = FALSE,     # Hide gene names (4,381 names will overlap and look messy)
  color = colorRampPalette(c("blue", "white", "red"))(50), # Blue=Low, White=0, Red=High
  main = "Hierarchical Clustering of Yeast Cell Cycle Genes"
)

# The "Spaghetti" Line Plot (Individual Trajectories)
# We pick just 5 specific genes and draw a line for each one over time.
# This lets us see the exact wave patterns of individual genes.
subset_5_genes <- unique(spellman_long$time)[1:5]
spellman_long %>%
  filter(time %in% subset_5_genes) %>%
  tidyplot(x = time_point, y = expression, color = time) %>%
  add_line(linewidth = 1) %>%
  add_data_points(size = 2) %>% # Adds dots at the actual measured time points
  adjust_size(width = 160, height = 80)

