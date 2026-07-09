# ==============================================================================
# SET THEORY COMPARISON: CORE vs. BRIDGE GENES (ROBUST VERSION)
# ==============================================================================

library(ggVennDiagram)
library(ggplot2)

# Load the lists directly from the CSV files we saved earlier
# (The $Gene_ID pulls just the column of names, converting it into a clean list)
super_core_253 <- read.csv("super_core_phase_5.csv")$Gene_ID
bridge_genes <- read.csv("umap_gmm_transient_bridges.csv")$Gene_ID

# Mathematical Set Operations
overlapping_genes <- intersect(super_core_253, bridge_genes)
pure_core <- setdiff(super_core_253, bridge_genes)
pure_bridge <- setdiff(bridge_genes, super_core_253)

# Print the formal comparison report
cat("\n--- Topological Gene Comparison ---\n")
cat("Total Phase 5 Super-Core Genes :", length(super_core_253), "\n")
cat("Total GMM Bridge Genes         :", length(bridge_genes), "\n")
cat("-----------------------------------\n")
cat("Overlap (Caught in BOTH lists) :", length(overlapping_genes), "genes\n")
cat("Pure Core (Only in Phase 5)    :", length(pure_core), "genes\n")
cat("Pure Bridge (Only in GMM borders):", length(pure_bridge), "genes\n")

# Visualizing the Comparison
gene_lists <- list(
  "Phase 5 Super-Core" = super_core_253,
  "Transient Bridges" = bridge_genes
)

# Plotting the Venn Diagram
ggVennDiagram(gene_lists) +
  scale_fill_gradient(low = "#F4FAFE", high = "#2B8CBE") +
  theme_void() +
  labs(
    title = "Topological Intersection: Core vs. Bridge Genes",
    subtitle = "Evaluating overlap between stable consensus centers and probabilistic boundaries"
  )
