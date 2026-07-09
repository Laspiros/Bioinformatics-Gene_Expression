# 1. Install the Bioconductor download manager from the standard store
 install.packages("BiocManager")

# 2. Use the BiocManager to download the elite biology packages
BiocManager::install(c("clusterProfiler", "org.Sc.sgd.db", "enrichplot"))

# ==============================================================================
# EXTERNAL BIOLOGICAL VALIDATION: GENE ONTOLOGY (GO) ENRICHMENT
# ==============================================================================

# 1. Install Bioconductor tools (Run these two lines ONLY if you don't have them yet)
# install.packages("BiocManager")
# BiocManager::install(c("clusterProfiler", "org.Sc.sgd.db", "enrichplot"))

library(clusterProfiler)
library(org.Sc.sgd.db)   # The official Yeast (Saccharomyces cerevisiae) database
library(ggplot2)
library(enrichplot)

cat("Running GO Enrichment Analysis. Querying the global yeast database...\n")

# 2. Load the elite genes you extracted (e.g., Phase 5)
phase_5_genes <- read.csv("super_core_phase_5.csv")$Gene_ID

# 3. Execute the GO Enrichment Algorithm
# We are testing for "BP" (Biological Process) to find out what these genes DO
go_results_phase5 <- enrichGO(
  gene          = phase_5_genes,
  OrgDb         = org.Sc.sgd.db,
  keyType       = "ORF",       # Standard Yeast Gene ID format (e.g., YAL001C)
  ont           = "BP",        # BP = Biological Process
  pAdjustMethod = "BH",        # Benjamini-Hochberg statistical correction
  pvalueCutoff  = 0.01,        # Only keep highly significant pathways
  qvalueCutoff  = 0.05
)

# 4. Generate the Industry-Standard "Dotplot"
dotplot(go_results_phase5, showCategory = 10) +
  labs(
    title = "Functional Validation: Phase 5 Super-Core",
    subtitle = "Gene Ontology (GO) Enrichment of the transient checkpoint vector",
    x = "Gene Ratio (Enrichment Strength)",
    y = "Biological Pathway"
  ) +
  theme_minimal()

