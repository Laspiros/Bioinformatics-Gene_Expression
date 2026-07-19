# Bioinformatics-Gene_Expression
Unsupervised ML pipeline isolating elite Phase 5 checkpoint genes in yeast.

# Unsupervised Machine Learning Pipeline for S. cerevisiae Cell Cycle

## Project Overview
This project maps the transcriptomic landscape of the Saccharomyces cerevisiae cell division cycle using a longitudinal dataset of 4,381 genes. By combining linear noise filtration (PCA) with a tri-model hard ensemble (K-Means, HC, PAM), the workflow isolates an absolute consensus framework of 2,953 Super-Core genes.

## Key Discoveries
* **The Phase 5 Checkpoint:** Intersecting the k=4 models revealed a hidden fifth sub-population (n=253).
* **Biological Function:** GO Enrichment Analysis proved Phase 5 is a rapid secretory checkpoint for ER membrane targeting (p = 0.0065) driven by 7 elite leader genes (including SGT2, GET4, and SEC63).

## Tech Stack
* **Language:** R
* **Libraries:** `cluster`, `mclust`, `ggplot2`, `org.Sc.sgd.db`, `clusterProfiler`
* **Techniques:** PCA, K-Means, Hierarchical Clustering, PAM, Gaussian Mixture Models (GMM)

## Files
* assets
* data
* scripts
* README
