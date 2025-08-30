# absolute path is safest; adjust folders to match your system
f <- "/Users/ben/Documents/UCRTP/Dutta Drive/output/A3SS.MATS.JCEC.HTMP.txt"

file.exists(f)            # should return TRUE
library(dplyr)

# 1. Read rMATS output (tab-separated)
a3   <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 2. Extract sample-level PSI columns (they start with IncLevel)
psi_cols <- grep("^IncLevel", names(a3))
psi_raw  <- a3[, psi_cols]

# 3. rMATS writes commas between biological replicates; split & average
split_and_mean <- function(x) {
  sapply(strsplit(x, ","), function(v) mean(as.numeric(v), na.rm = TRUE))
}

psi_matrix <- apply(psi_raw, 2, split_and_mean)              # matrix (events × samples)
rownames(psi_matrix) <- paste(a3$GeneID, a3$geneSymbol, sep = "_")
event_vars <- apply(psi_matrix, 1, var, na.rm = TRUE)
top_events <- order(event_vars, decreasing = TRUE)[1:50]
psi_top    <- psi_matrix[top_events, ]
colnames(psi_top)
sample_info <- data.frame(Group = factor(c(rep("Control", 3),
                                           rep("Knockdown", 2))))
# --- check the matrix first ---
colnames(psi_top)            # should print 6 names
ncol(psi_top)            # how many columns?  should be 6
length(colnames(psi_top))# should also be 6 and equal the value above

# If these are not both 6, stop here and figure out why.


# --- create a matching annotation data-frame ---
sample_info <- data.frame(
  Group = factor(c(rep("Control", 3),   # first 3 columns
                   rep("Knockdown", 3)) # last 3 columns
  )
)
rownames(sample_info) <- colnames(psi_top)   # lengths now match (6 = 6)

library(pheatmap)

pheatmap(
  psi_top,                          # 50×6 (or whatever) PSI matrix
  annotation_col = sample_info,     # 6-row annotation
  cluster_rows   = TRUE,
  cluster_cols   = TRUE,
  color          = colorRampPalette(c("royalblue", "white", "firebrick3"))(100),
  main           = "A3SS PSI Heatmap",
  fontsize_row   = 7,
  fontsize_col   = 10
)

