# ==============================================================================
# EXTRACTING THE 7 ELITE PATHWAY LEADER GENES
# ==============================================================================

# Convert the formal GO object into a standard R data frame
go_df <- as.data.frame(go_results_phase5)

# Filter the row corresponding exactly to the ER targeting pathway
er_pathway_row <- go_df[go_df$Description == "post-translational protein targeting to endoplasmic reticulum membrane", ]

# Extract the raw, slash-separated string of gene IDs
raw_gene_string <- er_pathway_row$geneID

# Clean and split the string into a pristine vector of individual gene names
checkpoint_leader_genes <- unlist(strsplit(raw_gene_string, "/"))

# --- Print the 7 genes directly to the console ---
cat("\n--- The 7 Elite Phase 5 Checkpoint Leader Genes ---\n")
print(checkpoint_leader_genes)

# Export these specific 7 genes to a dedicated CSV file for downstream biology
write.csv(data.frame(Leader_Gene_ID = checkpoint_leader_genes), 
          "phase_5_checkpoint_leaders.csv", 
          row.names = FALSE)

cat("\nSuccess: Saved 'phase_5_checkpoint_leaders.csv' to your directory!\n")