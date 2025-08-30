# Install packages if not already installed
install.packages("BiocManager")

BiocManager::install("maser")   # then run:  library(maser)
install.packages(c("ggplot2", "dplyr"))
# reinstall and then load
install.packages("ggplot2")   # or update.packages("ggplot2")
library(ggplot2)


# Load packages
library(maser)
library(ggplot2)
library(dplyr)

# Set path to your rMATS output directory
path_to_rmats <- "/Users/ben/Documents/UCRTP/Dutta Drive/output"  # note the **t** in output

# Create maser object with your condition labels
# Replace with your actual condition names
events <- maser(path_to_rmats,
                cond_labels = c("Control", "CARF_Knockdown"),
                ftype = "JCEC")        # or "JCEC"

events  

# Filter by coverage to remove low-confidence events
events_filtered <- filterByCoverage(events, avg_reads = 10)

# Create volcano plot for Skipped Exons (SE)
volcano_plot <- volcano(events_filtered, 
                        type = "SE",  # Can be "SE", "MXE", "A3SS", "A5SS", "RI"
                        fdr = 0.001,   # FDR threshold
                        deltaPSI = 0.2)  # PSI difference threshold

# Display the plot
print(volcano_plot)

# Save the plot
ggsave("volcano_plot_SE.png", volcano_plot, 
       width = 8, height = 6, dpi = 300)

