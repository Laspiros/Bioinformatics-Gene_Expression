install.packages("tidyplots")
install.packages("tidyverse")
install.packages("cluster")
install.packages("pheatmap")
install.packages("zoo")
install.packages("fpc")
install.packages("dbscan")


# Load the libraries into your script
library(tidyplots)
library(tidyverse)
library(cluster)
library(pheatmap)
library(zoo)
library(fpc)
library(ggforce)
library(mclust)
library(ggplot2)


# ===========================================================================
#  Data Description
# ===========================================================================

# Read the dataset
spellman_data <- read.csv("Spellman.csv")
# Explore the dataset
dim(spellman_data)
str(spellman_data)
head(spellman_data)

# ===========================================================================
#  Data Cleaning - EDA (Exploratory Data Analysis)
# ===========================================================================

# Handle Missing Values (NAs)
# Look at how many missing values we have before cleaning
total_nas <- sum(is.na(spellman_data))
cat("Total missing values in original data:", total_nas, "\n")

# Set Row Names
# Convert column to row names
spellman_data2 <- spellman_data
rownames(spellman_data2) <- spellman_data$time
spellman_data2$time <- NULL # Remove the old column

# ===========================================================================
#  Check data statistics
# ===========================================================================

# Convert the data frame to a pure numeric matrix
spellman_matrix <- as.matrix(spellman_data2)

# Calculate the standard deviation
global_sd <- sd(spellman_matrix)

# Print it out to verify
print(global_sd)

# Calculate the rest of descriptive statistics 
global_mean <- mean(spellman_matrix)
global_min  <- min(spellman_matrix)
global_max  <- max(spellman_matrix)

# Outlier Detection Rules
upper_bound <- global_mean + (3 * global_sd)
lower_bound <- global_mean - (3 * global_sd)

# Count how many data points are outside these bounds
outliers_count <- sum(spellman_matrix > upper_bound | spellman_matrix < lower_bound)
total_data_points <- length(spellman_matrix)
outlier_percentage <- (outliers_count / total_data_points) * 100

# Print the report
cat("========== DATA STATISTICS REPORT ==========\n")
cat("Absolute Minimum Expression: ", global_min, "\n")
cat("Absolute Maximum Expression: ", global_max, "\n")
cat("Overall Average (Mean):      ", global_mean, "\n")
cat("Upper Outlier Threshold (+3 SD): ", upper_bound, "\n")
cat("Lower Outlier Threshold (-3 SD): ", lower_bound, "\n")
cat("Percentage of Data as Outliers:  ", round(outlier_percentage, 2), "%\n")

