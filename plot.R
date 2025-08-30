install.packages("BiocManager")
BiocManager::install("maser")
install.packages(c("ggplot2", "dplyr"))

# Load packages
library(maser)
library(ggplot2)
library(dplyr)


# Set path to your rMATS output directory
path_to_rmats <- path_to_rmats <- "/Users/ben/Documents/UCRTP/Dutta Drive/output"

# Create maser object with your condition labels
# Replace with your actual condition names
events <- maser(path_to_rmats, 
                cond_labels = c("Control", "CARF_Knockdown"),
                ftype = "JCEC")  # or "JCEC" depending on your files

# Filter by coverage to remove low-confidence events
events_filtered <- filterByCoverage(events, avg_reads = 10)


# Get top significant events
top_events <- topEvents(events_filtered, 
                        type = "SE", 
                        fdr = 0.05, 
                        deltaPSI = 0.1)

# Extract PSI values for these events
psi_matrix <- PSI(top_events, type = "SE")

# Take top 50 events for visualization
psi_top50 <- psi_matrix[1:50, ]
