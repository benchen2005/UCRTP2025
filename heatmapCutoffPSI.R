###############################################################################
# USER SETTINGS                                                               #
###############################################################################
file_path  <- "/Users/ben/Documents/Filtered SE.tsv"                 # TSV path
psi_cols   <- c("Ctrl1","Ctrl2","Ctrl3","KD1","KD2","KD3")           # PSI cols
psi_cutoff <- 0.20                                                   # |mean PSI|
out_dir    <- "plots"
out_png    <- file.path(out_dir, "SE_PSI_heatmap_absPSI>=0.2.png")

###############################################################################
# PACKAGES (auto-install if missing)                                          #
###############################################################################
pkgs <- c("dplyr","readr","pheatmap")
need <- setdiff(pkgs, rownames(installed.packages()))
if (length(need)) install.packages(need)
lapply(pkgs, library, character.only = TRUE)

###############################################################################
# READ TSV & BUILD PSI MATRIX                                                 #
###############################################################################
stopifnot(file.exists(file_path))
se <- read.delim(file_path, check.names = FALSE, stringsAsFactors = FALSE)

missing <- setdiff(psi_cols, names(se))
if (length(missing))
  stop("Missing PSI columns: ", paste(missing, collapse = ", "))

psi_df <- se %>%
  mutate(event_id = paste0(geneSymbol, "_", row_number())) %>%     # unique row id
  select(event_id, all_of(psi_cols))

rownames(psi_df) <- psi_df$event_id
psi_df$event_id  <- NULL
psi_mat <- as.matrix(psi_df); mode(psi_mat) <- "numeric"

###############################################################################
# FILTER: keep events with |mean PSI| ≥ 0.20                                  #
###############################################################################
mean_psi <- rowMeans(psi_mat, na.rm = TRUE)
keep     <- abs(mean_psi) >= psi_cutoff
psi_mat  <- psi_mat[keep, ]

cat("Retained", nrow(psi_mat), "events with |mean PSI| ≥", psi_cutoff, "\n")

###############################################################################
# SAMPLE ANNOTATION                                                           #
###############################################################################
groups <- ifelse(grepl("^Ctrl", colnames(psi_mat)), "Control", "KnockDown")
anno   <- data.frame(Group = factor(groups, levels = c("Control","KnockDown")))
rownames(anno) <- colnames(psi_mat)

###############################################################################
# ENSURE OUTPUT DIR                                                           #
###############################################################################
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

###############################################################################
# DRAW & SAVE HEAT-MAP                                                        #
###############################################################################
png(out_png, width = 2400, height = 1800, res = 300)
pheatmap::pheatmap(
  psi_mat,
  annotation_col = anno,
  show_rownames  = FALSE,
  cluster_rows   = TRUE,
  cluster_cols   = TRUE,
  color          = colorRampPalette(c("royalblue","white","firebrick3"))(100),
  main           = sprintf("PSI heat-map | retained %d events (|mean PSI| ≥ %.2f)",
                           nrow(psi_mat), psi_cutoff)
)
dev.off()

cat("✔ Heat-map saved to:", normalizePath(out_png), "\n")

