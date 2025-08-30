###############################################################################
# 0 — USER SETTINGS                                                           #
###############################################################################
file_path  <- "/Users/ben/Documents/UCRTP/Filtered Events/noFilterFDR.SE.tsv"
psi_cols   <- c("Ctrl1","Ctrl2","Ctrl3","KD1","KD2","KD3")  # PSI columns (must exist)
out_file   <- "heatmap_SRSF_family.png"

# SRSF family + common aliases → map all aliases to canonical gene symbol
srsf_aliases <- list(
  SRSF1 = c("SRSF1","SF2","ASF","ASF/SF2"),
  SRSF2 = c("SRSF2","SC35","SRp30b"),
  SRSF3 = "SRSF3",
  SRSF4 = "SRSF4",
  SRSF5 = "SRSF5",
  SRSF6 = "SRSF6",
  SRSF7 = "SRSF7",
  SRSF10 = "SRSF10"
)

###############################################################################
# 1 — PACKAGES                                                                #
###############################################################################
pkgs <- c("dplyr","readr","pheatmap")
new  <- setdiff(pkgs, rownames(installed.packages()))
if (length(new)) install.packages(new)
invisible(lapply(pkgs, library, character.only = TRUE))

# helper: safe numeric parse (handles characters like "0.45" or "0.45%")
numify <- function(x) if (is.character(x)) readr::parse_number(x) else as.numeric(x)

###############################################################################
# 2 — LOAD DATA & CHECKS                                                      #
###############################################################################
stopifnot(file.exists(file_path))
se <- utils::read.delim(file_path, check.names = FALSE, stringsAsFactors = FALSE)

if (!"geneSymbol" %in% names(se)) {
  stop("Column 'geneSymbol' not found in the TSV: ", basename(file_path))
}

missing <- setdiff(psi_cols, names(se))
if (length(missing)) {
  stop("Missing PSI columns: ", paste(missing, collapse = ", "),
       "\nPresent columns are: ", paste(names(se), collapse = ", "))
}

# build alias → canonical map (case-insensitive match)
alias_map <- unlist(srsf_aliases, use.names = FALSE)
canon_vec <- rep(names(srsf_aliases), lengths(srsf_aliases))
names(canon_vec) <- tolower(alias_map)  # key = alias (lowercase), value = canonical

# annotate canonical SRSF name if gene matches any alias
se$._gene_lower <- tolower(se$geneSymbol)
se$SRSF_canonical <- canon_vec[se$._gene_lower]
srsf_se <- se[!is.na(se$SRSF_canonical), , drop = FALSE]

if (nrow(srsf_se) == 0L) {
  stop("No SRSF family hits found.\nChecked aliases: ",
       paste(unique(alias_map), collapse = ", "))
}

###############################################################################
# 3 — PSI MATRIX                                                              #
###############################################################################
# keep PSI columns and create unique event IDs labeled by canonical gene
psi_df <- srsf_se |>
  dplyr::mutate(event_id = paste0(SRSF_canonical, "_", dplyr::row_number())) |>
  dplyr::select(event_id, dplyr::all_of(psi_cols))

# rownames as unique IDs, then drop the id column
rownames(psi_df) <- psi_df$event_id
psi_df$event_id  <- NULL

# coerce PSI columns to numeric
psi_df[] <- lapply(psi_df, numify)

psi_mat <- as.matrix(psi_df)  # events × samples

if (anyNA(psi_mat)) {
  message("Note: NA values present in PSI matrix (will appear as gaps in the heatmap).")
}

###############################################################################
# 4 — SAMPLE ANNOTATION                                                       #
###############################################################################
sample_groups <- ifelse(grepl("^Ctrl", colnames(psi_mat), ignore.case = FALSE),
                        "Control", "KnockDown")
anno <- data.frame(Group = factor(sample_groups, levels = c("Control","KnockDown")))
rownames(anno) <- colnames(psi_mat)

###############################################################################
# 5 — DRAW + SAVE HEATMAP (png device to avoid ComplexHeatmap filename issue) #
###############################################################################
# choose a diverging palette for PSI (0..1); if you prefer ΔPSI (−1..1), adjust upstream
pal <- colorRampPalette(c("royalblue","white","firebrick3"))(100)

# on-screen preview
pheatmap::pheatmap(
  psi_mat,
  annotation_col = anno,
  cluster_rows   = TRUE,
  cluster_cols   = TRUE,
  color          = pal,
  main           = sprintf("PSI heatmap — SRSF family (%d events)", nrow(psi_mat)),
  fontsize_row   = 11,
  fontsize_col   = 10,
  legend         = TRUE,
  border_color   = NA
)

# save via png() to be robust across environments
dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
png(filename = out_file, width = 8*300, height = 10*300, res = 300)
pheatmap::pheatmap(
  psi_mat,
  annotation_col = anno,
  cluster_rows   = TRUE,
  cluster_cols   = TRUE,
  color          = pal,
  main           = sprintf("PSI heatmap — SRSF family (%d events)", nrow(psi_mat)),
  fontsize_row   = 14,
  fontsize_col   = 10,
  legend         = TRUE,
  border_color   = NA
)
dev.off()

cat("✔ Heatmap saved to", normalizePath(out_file), "\n")

