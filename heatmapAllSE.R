###############################################################################
# USER SETTINGS                                                               #
###############################################################################
file_path  <- "/Users/ben/Documents/Filtered SE.tsv"               # TSV file
psi_cols   <- c("Ctrl1","Ctrl2","Ctrl3","KD1","KD2","KD3")         # PSI columns
row_mode   <- "ALL"          # "ALL" | "TOP_VAR" | "PSI_CUTOFF"
top_n      <- 1000           # used only if row_mode == "TOP_VAR"
psi_cutoff <- 0.80           # used only if row_mode == "PSI_CUTOFF"
out_dir    <- "plots"        # output folder

out_png <- file.path(
  out_dir,
  switch(row_mode,
         ALL        = "SE_PSI_heatmap_ALL.png",
         TOP_VAR    = paste0("SE_PSI_heatmap_top", top_n, ".png"),
         PSI_CUTOFF = paste0("SE_PSI_heatmap_absPSI≥", psi_cutoff, ".png"))
)

###############################################################################
# PACKAGES (auto-install if missing)                                          #
###############################################################################
pkgs <- c("dplyr","tidyr","readr","pheatmap")
need <- setdiff(pkgs, rownames(installed.packages()))
if (length(need)) install.packages(need)
lapply(pkgs, library, character.only = TRUE)

###############################################################################
# READ TSV & VERIFY PSI COLUMNS                                               #
###############################################################################
stopifnot(file.exists(file_path))
se <- read.delim(file_path, check.names = FALSE, stringsAsFactors = FALSE)

miss <- setdiff(psi_cols, names(se))
if (length(miss))
  stop("Missing PSI columns: ", paste(miss, collapse = ", "))

###############################################################################
# BUILD NUMERIC PSI MATRIX (events × samples)                                 #
###############################################################################
psi_df <- se %>%
  mutate(event_id = paste0(geneSymbol, "_", row_number())) %>%  # unique ID
  select(event_id, all_of(psi_cols))

rownames(psi_df) <- psi_df$event_id
psi_df$event_id  <- NULL
psi_mat <- as.matrix(psi_df); mode(psi_mat) <- "numeric"

###############################################################################
# ROW-SELECTION LOGIC                                                         #
###############################################################################
if (row_mode == "TOP_VAR") {
  vars <- apply(psi_mat, 1, var, na.rm = TRUE)
  keep <- order(vars, decreasing = TRUE)[1:min(top_n, nrow(psi_mat))]
  psi_mat <- psi_mat[keep, ]
  message("Keeping top ", nrow(psi_mat), " events by variance")
} else if (row_mode == "PSI_CUTOFF") {
  mean_psi <- rowMeans(psi_mat, na.rm = TRUE)
  keep     <- abs(mean_psi) >= psi_cutoff
  psi_mat  <- psi_mat[keep, ]
  message("Keeping ", nrow(psi_mat),
          " events with |mean PSI| ≥ ", psi_cutoff)
} else {
  message("Keeping ALL ", nrow(psi_mat), " events")
}

###############################################################################
# SAMPLE ANNOTATION (Control vs KD)                                           #
###############################################################################
sample_groups <- ifelse(grepl("^Ctrl", colnames(psi_mat)), "Control","KnockDown")
anno <- data.frame(Group = factor(sample_groups, levels = c("Control","KnockDown")))
rownames(anno) <- colnames(psi_mat)

###############################################################################
# ENSURE OUTPUT DIRECTORY                                                     #
###############################################################################
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

###############################################################################
# DRAW HEAT-MAP ON SCREEN                                                     #
###############################################################################
pheatmap::pheatmap(
  psi_mat,
  annotation_col = anno,
  show_rownames  = FALSE,
  cluster_rows   = TRUE,
  cluster_cols   = TRUE,
  color          = colorRampPalette(c("royalblue","white","firebrick3"))(100),
  main           = sprintf("PSI heat-map (%s rows, %s samples)",
                           nrow(psi_mat), ncol(psi_mat))
)

###############################################################################
# SAVE PNG (works even if ComplexHeatmap is loaded)                           #
###############################################################################
png(out_png, width = 2400, height = 1800, res = 300)
pheatmap::pheatmap(
  psi_mat,
  annotation_col = anno,
  show_rownames  = FALSE,
  cluster_rows   = TRUE,
  cluster_cols   = TRUE,
  color          = colorRampPalette(c("royalblue","white","firebrick3"))(100),
  main           = sprintf("PSI heat-map (%s rows, %s samples)",
                           nrow(psi_mat), ncol(psi_mat))
)
dev.off()

cat("✔ Heat-map saved to", normalizePath(out_png), "\n")

